module mynn

! References: "NN09" Nakanishi, M. and Niino, 2009: Development of an improved turbulence closure model
!                            for the atmospheric boundary layer. J. Meteor. Soc. Japan, 87, 895-912

use MAPL_ConstantsMod, only: MAPL_GRAV, MAPL_CP, MAPL_ALHL, MAPL_RGAS, MAPL_RVAP, MAPL_P00, MAPL_VIREPS

implicit none

real, parameter :: onethird  = 1./3.
real, parameter :: twothirds = 2./3.

real, parameter :: gocp  = MAPL_GRAV/MAPL_CP
real, parameter :: lvocp = MAPL_ALHL/MAPL_CP
real, parameter :: kappa = MAPL_RGAS/MAPL_CP
real, parameter :: ep2   = MAPL_RVAP/MAPL_RGAS - 1.

! MYNN constants 
real, parameter :: Pr     = 0.74
real, parameter :: gamma1 = 0.235
real, parameter :: B1     = 24.
real, parameter :: B2     = 15.
real, parameter :: C2     = 0.729
real, parameter :: C4     = 0.
real, parameter :: C5     = 0.2

double precision, parameter :: C3 = 0.34d0

real, parameter :: A1 = B1*( 1. - 3.*gamma1 )/6.

double precision, parameter :: C1 = gamma1 - 1./(3.*A1*B1**onethird)

real, parameter :: A2     = A1*( gamma1 - C1 )/(gamma1*Pr)
real, parameter :: gamma2 = B2/B1*( 1. - C3 ) + 2.*A1/B1*( 3. - 2.*C2 )

double precision, parameter :: e1c = 3.*A2*B2*( 1. - C3 )
double precision, parameter :: e2c = 9.*A1*A2*( 1. - C2 )
double precision, parameter :: e3c = 9.*A2**2.*( 1. - C2 )*( 1. - C5 )
double precision, parameter :: e4c = 12.*A1*A2*( 1. - C2 )
double precision, parameter :: e5c = 6.*A1**2.

real, parameter :: eMc = 3.*A1*( 1. - C3 )
real, parameter :: eHc = 3.*A2*( 1. - C3 )

real, parameter :: F1 = B1*( gamma1 - C1 ) + 2.*A1*( 3. - 2.*C2 ) + 3.*A2*( 1. - C2 )*( 1. - C5 ) ! NN09 (A6)
real, parameter :: F2 = B1*( gamma1 + gamma2 ) - 3.*A1*( 1. - C2 )                                ! NN09 (A7)

double precision, parameter :: Rfc = gamma1/( gamma1 + gamma2 )                                    ! NN09 (A10)
double precision, parameter :: Rf1 = B1*( gamma1 - C1 )/F1                                         ! NN09 (A8)
double precision, parameter :: Rf2 = B1*gamma1/F2                                                  ! NN09 (A9)

real, parameter :: SMc = A1/A2*F1/F2 
real, parameter :: SHc = 3.*A2*( gamma1 + gamma2 )

double precision, parameter :: Ri1 = 0.5/SMc             ! NN09 (A12)
double precision, parameter :: Ri2 = Rf1*SMc             ! NN09 (A13)
double precision, parameter :: Ri3 = 4.*Rf2*SMc - 2.*Ri2 ! NN09 (A14)
double precision, parameter :: Ri4 = Ri2**2.

real, parameter :: c_tke_surf = B1**(2./3.)
real, parameter :: c_th2_surf = B2/B1**(1./3.)


! Test
logical :: initialized_mynn = .false.

contains

!
! run_mynn
!
subroutine run_mynn(IM, JM, LM, &                                                           ! in
                    DEBUG_FLAG, DOMF, MYNN_LEVEL, CONSISTENT_TYPE, &                        ! in
                    local_flag, WQL_TYPE, WRF_CG_FLAG, &                                    ! in
                    alpha1, alpha2, alpha3, alpha4, &                                       ! in (length scale parameters)
                    th00, ice_ramp, ple, pl, rhoe, zle, zlo, &                              ! in
                    u, v, T, qv, ql, qi, thl, qt, thv, &                                    ! in (mean atmospheric state)
                    u_star, H_surf, E_surf, &                                               ! in (surface state)
                    whl_mf, wqt_mf, wthv_mf, au, Mu, wu, E, D, wdet, &                      ! in (updraft state)
                    aci, Ai_moist, Bi_moist, &                                              ! in
                    tke, hl2, qt2, hlqt, &                                                  ! inout (turbulent state)
                    ws_explicit, wqv_explicit, wql_explicit, &                              ! inout (counter-gradient fluxes)
                    KM, KH, K_tke, itau, qdiv_out, SM25, SH25, L_out, &                     ! out
                    Beta_hl, Beta_qt, EM, EH, &                                             ! out
                    LS_out, LB_out, LT_out, &                                               ! out
                    tket_M, tket_B, tket_T_mf, hl2t_M, qt2t_M, hlqtt_M, &                   ! out
                    tket_T_mf1, tket_T_mf2, tket_T_mf3, tket_T_mf4, &                       ! out
                    tke_surf, hl2_surf, qt2_surf, hlqt_surf)                                ! out

  use MAPL_ConstantsMod, only: MAPL_KARMAN
  use MAPL_SatVaporMod, only: MAPL_EQsat

  integer, intent(in)                        :: IM, JM, LM, MYNN_LEVEL, CONSISTENT_TYPE, &
                                                WQL_TYPE, WRF_CG_FLAG, DEBUG_FLAG, local_flag
  real, intent(in)                           :: th00, ice_ramp, DOMF
  double precision, intent(in)               :: alpha1, alpha2, alpha3, alpha4
  real, dimension(IM,JM), intent(in)         :: u_star, H_surf, E_surf
  real, dimension(IM,JM,LM), intent(in)      :: zlo, u, v, T, qv, ql, qi, thv, thl, qt, E, D, wdet, pl
  real, dimension(IM,JM,0:LM), intent(in)    :: ple, zle, rhoe, whl_mf, wqt_mf, wthv_mf, &
                                                au, Mu, wu, Ai_moist, Bi_moist, aci
  real, dimension(IM,JM,0:LM), intent(inout) :: tke, hl2, qt2, hlqt, ws_explicit, wqv_explicit, wql_explicit
  real, dimension(IM,JM), intent(out)        :: tke_surf, hl2_surf, qt2_surf, hlqt_surf, LT_out
  real, dimension(IM,JM,0:LM), intent(out)   :: KM, KH, itau, qdiv_out, SM25, SH25, EM, EH, L_out, Beta_hl, Beta_qt, &
                                                tket_M, tket_B, tket_T_mf, hl2t_M, qt2t_M, hlqtt_M, &
                                                tket_T_mf1, tket_T_mf2, tket_T_mf3, tket_T_mf4, &
                                                LS_out, LB_out
  real, dimension(IM,JM,LM), intent(out)     :: K_tke

  integer :: i, j, k, kp1, km1

  double precision :: GH, GM, q2, Phi1, Phi2, Phi3, Phi4, Phi5, &
                      D_25, qdiv, qdiv2, L2, L2GM, L2GH, Lq, &
                      tau_test, w2_test, wb_test_1, wb_test_2 ! for debugging

  real :: dhldz, dqtdz, dqldz, idzlo, iexner, &
          whl, wqt, th_star, q_star, wrk1, &
          whl_explicit, wqt_explicit, wb_explicit, &
          hl2_25, qt2_25, hlqt_25, &
          goth00, Cw_low, Cw_high, dzle, &
          hlthv_25, qtthv_25, thv2_25 ! for debugging2

  real, dimension(IM,JM)      :: wb_surf, LMO
  real, dimension(IM,JM,LM)   :: hl
  real, dimension(IM,JM,0:LM) :: S2, N2, zeta, A, B

  double precision, dimension(IM,JM)      :: w_star, LT
  double precision, dimension(IM,JM,0:LM) :: q, L, LS, LB

  goth00 = MAPL_GRAV/th00

  ! Compute conserved thermodynamic properties
  do k = 1,LM
     do j = 1,JM
     do i = 1,IM
        hl(i,j,k) = T(i,j,k) + gocp*zlo(i,j,k) - lvocp*ql(i,j,k)
     end do
     end do
  end do

  ! Compute shear and buoyancy frequencies
  do k = 1,LM-1

     kp1 = k + 1
     do j = 1,JM
     do i = 1,IM
        q(i,j,k) = sqrt( max( 1.d-10, real(2.*tke(i,j,k),8) ) )

        idzlo = 1./( zlo(i,j,k) - zlo(i,j,kp1) )
        dhldz = ( hl(i,j,k) - hl(i,j,kp1) )*idzlo
        dqtdz = ( qt(i,j,k) - qt(i,j,kp1) )*idzlo        
        
        iexner = (MAPL_P00/ple(i,j,k))**kappa
        wrk1   = lvocp*iexner - ( 1. + ep2 )*th00

        Beta_hl(i,j,k) = iexner*( 1. - wrk1*aci(i,j,k)*Bi_moist(i,j,k) )
        Beta_qt(i,j,k) = ep2*th00 + wrk1*aci(i,j,k)*Ai_moist(i,j,k)

        N2(i,j,k) = goth00*( Beta_hl(i,j,k)*dhldz + Beta_qt(i,j,k)*dqtdz )
        S2(i,j,k) = (( u(i,j,k) - u(i,j,kp1) )*idzlo)**2. + (( v(i,j,k) - v(i,j,kp1) )*idzlo)**2.
     end do
     end do
  end do

  ! Compute some surface and TOA quantities
  do j = 1,JM
  do i = 1,IM
     wb_surf(i,j) = goth00*( H_surf(i,j)/(MAPL_CP*rhoe(i,j,LM)) + ep2*th00*E_surf(i,j)/rhoe(i,j,LM) )
     LMO(i,j)     = -u_star(i,j)**3./(MAPL_KARMAN*wb_surf(i,j))

     KM(i,j,0)  = 0.
     KM(i,j,LM) = 0.

     KH(i,j,0)  = 0.
     KH(i,j,LM) = 0.

     q(i,j,0)  = 0.d0
     q(i,j,LM) = 0.d0

     N2(i,j,0)  = 0.
     N2(i,j,LM) = 0.

     S2(i,j,0)  = 0.
     S2(i,j,LM) = 0.

     tket_B(i,j,LM) = wb_surf(i,j)

     th_star = H_surf(i,j)/(MAPL_CP*rhoe(i,j,LM)*u_star(i,j))
     q_star  = E_surf(i,j)/(rhoe(i,j,LM)*u_star(i,j))

     tke_surf(i,j)  = 0.5*c_tke_surf*u_star(i,j)**2.
     hl2_surf(i,j)  = c_th2_surf*th_star**2.
     qt2_surf(i,j)  = c_th2_surf*q_star**2.
     hlqt_surf(i,j) = sqrt(hl2_surf(i,j)*qt2_surf(i,j)) ! temporary
  end do
  end do

  ! Test
  if ( .not. initialized_mynn ) then
     call initialize_mynn(IM, JM, LM, &
                          local_flag, alpha1, alpha2, alpha3, alpha4, &
                          th00, ice_ramp, hl, qt, tke, hl2, qt2, hlqt, q, &
                          zle, zlo, S2, N2, &
                          u_star, wb_surf, LMO, &
                          thv, ple, pl)

     initialized_mynn = .true.
  end if

  call mynn_length(IM, JM, LM, &                                    ! in
                   local_flag, alpha1, alpha2, alpha3, alpha4, &    ! in
                   th00, ice_ramp, wb_surf, zle, zlo, q, N2, LMO, & ! in
                   thv, ple, pl, &                                  ! in
                   L, LS, LB, LT, w_star)                           ! out

  do k = 1,LM-1

     km1 = k - 1
     kp1 = k + 1
     do j = 1,JM
     do i = 1,IM
        ! Compute things needed for discretization
        idzlo   = 1./( zlo(i,j,k) - zlo(i,j,kp1) )
        dhldz   = ( hl(i,j,k) - hl(i,j,kp1) )*idzlo
        dqtdz   = ( qt(i,j,k) - qt(i,j,kp1) )*idzlo
        dqldz   = ( ql(i,j,k) - ql(i,j,kp1) )*idzlo

        dzle = zle(i,j,km1) - zle(i,j,k)

        ! Diagnose MYNN level-2.5 state
        call mynn_l25(L(i,j,k), q(i,j,k), tke(i,j,k), dhldz, dqtdz, S2(i,j,k), N2(i,j,k), & ! in
                      SM25(i,j,k), SH25(i,j,k), &                                           ! out
                      Phi1, Phi2, Phi3, Phi4, Phi5, D_25, L2, GH, L2GH, q2, qdiv, qdiv2)    ! out (needed for MYNN Level-3)

        ! 
        qdiv_out(i,j,k) = qdiv

        !
        itau(i,j,k) = q(i,j,k)/L(i,j,k)
        Lq          = L(i,j,k)*real( sqrt(2.*tke(i,j,k)), 8 )

        ! Compute thermodyanamic (co-)variances from level-2.5 closure
        hl2_25  = qdiv*B2*L2*SH25(i,j,k)*dhldz**2.
        qt2_25  = qdiv*B2*L2*SH25(i,j,k)*dqtdz**2.
        hlqt_25 = qdiv*B2*L2*SH25(i,j,k)*dhldz*dqtdz

        !
        if ( MYNN_LEVEL == 2 ) then
           ! Update thermodynamic second-order moments
           hl2(i,j,k)  = hl2_25
           qt2(i,j,k)  = qt2_25
           hlqt(i,j,k) = hlqt_25
        else           
           ! 
           if ( MYNN_LEVEL == 4 ) then
              hl2(i,j,k)  = hl2_25
              hlqt(i,j,k) = hlqt_25
           end if

           ! Compute MYNN level-3 state
           call mynn_l3(th00, Phi1, Phi3, Phi4, Phi5, D_25, L2, GH, L2GH, qdiv, qdiv2, & ! in (MYNN level-2.5 state)
                        Phi2, q2, &                                                      ! inout          
                        EM(i,j,k), EH(i,j,k), Cw_low, Cw_high)                           ! out (MYNN level-3 state)
        end if

        ! Compute diffusivities and second-order tendencies
        call mynn_tendency(mynn_level, wrf_cg_flag, domf, consistent_type, &                                               ! in
                           th00, rhoe(i,j,k), Beta_hl(i,j,k), Beta_qt(i,j,k), hl2(i,j,k), qt2(i,j,k), hlqt(i,j,k), &       ! in
                           dzle, dhldz, dqtdz, dqldz, N2(i,j,k), S2(i,j,k), &                                              ! in
                           L2, Lq, SM25(i,j,k), SH25(i,j,k), EM(i,j,k), EH(i,j,k), hl2_25, qt2_25, hlqt_25, &              ! in
                           tke(i,j,k), tke(i,j,km1), Mu(i,j,k), Mu(i,j,km1), au(i,j,k), au(i,j,km1), &                     ! in (for tket_T_mf)
                           wu(i,j,k), wu(i,j,km1), E(i,j,k), D(i,j,k), wdet(i,j,k), &                                      ! in (for tket_T_mf)
                           whl_mf(i,j,k), wqt_mf(i,j,k), wthv_mf(i,j,k), &                                                 ! in
                           KM(i,j,k), KH(i,j,k), ws_explicit(i,j,k), wqv_explicit(i,j,k), wql_explicit(i,j,k), &           ! out
                           tket_M(i,j,k), tket_B(i,j,k), tket_t_mf(i,j,k), hl2t_M(i,j,k), qt2t_M(i,j,k), hlqtt_M(i,j,k), & ! out
                           tket_T_mf1(i,j,k), tket_T_mf2(i,j,k), tket_T_mf3(i,j,k), tket_T_mf4(i,j,k))                     ! out

        ! Debugging output
        if ( DEBUG_FLAG == 1 ) then
           tau_test = L(i,j,k)/q(i,j,k)
           w2_test  = onethird*q(i,j,k)**2. + 2.*A1*tau_test*( -KM(i,j,k)*S2(i,j,k) ) &
                                            + 4.*A1*( 1. - C2 )*tau_test*( -KH(i,j,k)*N2(i,j,k) + wb_explicit)
           if ( MYNN_LEVEL == 2 ) then
              hlthv_25 = Beta_hl(i,j,k)*hl2_25  + Beta_qt(i,j,k)*hlqt_25
              qtthv_25 = Beta_hl(i,j,k)*hlqt_25 + Beta_qt(i,j,k)*qt2_25
              thv2_25  = max( 0.d0, Beta_hl(i,j,k)*hlthv_25 + Beta_qt(i,j,k)*qtthv_25 )

              wb_test_1 = MAPL_CP*rhoe(i,j,k)*3.*A2*tau_test*( -w2_test*N2(i,j,k) + ( 1. - C3 )*goth00**2.*thv2_25 )
              wb_test_2 = MAPL_CP*rhoe(i,j,k)*( -KH(i,j,k)*N2(i,j,k) + wb_explicit )
           else
!              wb_test  = 3.*A2*tau_test*(-w2_test*N2(i,j,k) + ( 1. - C3 )*goth002*( thv2_25  + thv2_p ) )
           end if

           write(*,*)             &
                      tke(i,j,k), &
!                      wb_test_1, &
!                      wb_test_2, &
                      aci(i,j,k), &
                      real(qdiv, 4)
        end if
     end do
     end do
  end do

  ! Compute TKE diffusivity
  do k = 2,LM-1
     km1 = k - 1
     do j = 1,JM
     do i = 1,IM
        K_tke(i,j,k) = 1.5*( KM(i,j,k) + KM(i,j,km1) )
     end do
     end do
  end do
  K_tke(:,:,1)  = K_tke(:,:,2)
  K_tke(:,:,LM) = K_tke(:,:,LM-1)

  !
  L_out  = L
  LS_out = LS
  LB_out = LB
  LT_out = LT

end subroutine run_mynn

!
! mynn_tendency
!
subroutine mynn_tendency(mynn_level, wrf_cg_flag, domf, consistent_type, &      ! in
                         th00, rhoe, Beta_hl, Beta_qt, hl2, qt2, hlqt, &        ! in
                         dzle, dhldz, dqtdz, dqldz, N2, S2, &                   ! in
                         L2, Lq, SM25, SH, EM, EH, hl2_25, qt2_25, hlqt_25, &   ! in
                         tke, tke_up, Mu, Mu_up, au, au_up, &                   ! in
                         wu, wu_up, E, D, wdet, &                               ! in
                         whl_mf, wqt_mf, wthv_mf, &                             ! in
                         KM, KH, ws_explicit, wqv_explicit, wql_explicit, &     ! out
                         tket_M, tket_B, tket_t_mf, hl2t_M, qt2t_M, hlqtt_M, &  ! out
                         tket_T_mf1, tket_T_mf2, tket_T_mf3, tket_T_mf4)        ! out

  implicit none

  integer, intent(in)          :: mynn_level, consistent_type, wrf_cg_flag
  real, intent(in)             :: domf, th00, rhoe, Beta_hl, Beta_qt, dzle, dhldz, dqtdz, dqldz, N2, S2, &
                                  SM25, SH, EM, EH, &
                                  whl_mf, wqt_mf, wthv_mf, hl2, qt2, hlqt, hl2_25, qt2_25, hlqt_25, &
                                  tke, tke_up, Mu, Mu_up, au, au_up, wu, wu_up, E, D, wdet
  double precision, intent(in) :: L2, Lq
  real, intent(out)            :: KM, KH, ws_explicit, wqv_explicit, wql_explicit, &                                  
                                  tket_M, tket_B, tket_t_mf, hl2t_M, qt2t_M, hlqtt_M, &
                                  tket_T_mf1, tket_T_mf2, tket_T_mf3, tket_T_mf4

  real :: goth00, goth002, Cw_low, Cw_high, Cw_25, hlthv, qtthv, thv2, &
          hlthv_p, qtthv_p, thv2_p, SM, whl_explicit, wqt_explicit, wb_explicit, &
          whl, wqt, hlthv_25, qtthv_25, thv2_25, we, we_up, tkee, tkee_up, w2e_up

  goth00  = MAPL_GRAV/th00
  goth002 = goth002**2.

  if ( mynn_level >= 3 ) then
     ! Compute level-2.5 buoyancy (co-)variances
     hlthv_25 = Beta_hl*hl2_25  + Beta_qt*hlqt_25
     qtthv_25 = Beta_hl*hlqt_25 + Beta_qt*qt2_25
     thv2_25  = max( 0., Beta_hl*hlthv_25 + Beta_qt*qtthv_25 )

     ! Compute bouyancy (co-variances)
     hlthv = Beta_hl*hl2  + Beta_qt*hlqt
     qtthv = Beta_hl*hlqt + Beta_qt*qt2
     thv2  = max( 0.d0, Beta_hl*hlthv + Beta_qt*qtthv )

     ! Compute counter-gradient flux of conserved variables
     if ( wrf_cg_flag == 0 ) then
        hlthv_p = hlthv - hlthv_25
        qtthv_p = qtthv - qtthv_25
     else
        if ( hlthv_25 >= 0. ) then
           hlthv_p = max( 0., hlthv - hlthv_25 )
        else
           hlthv_p = min( 0., hlthv - hlthv_25 )
        end if
        
        if ( qtthv_25 >= 0. ) then
           qtthv_p = max( 0., qtthv - qtthv_25 )
        else
           qtthv_p = min( 0., qtthv - qtthv_25 )
        end if
     end if
     whl_explicit = Lq*EH*goth00*hlthv_p
     wqt_explicit = Lq*EH*goth00*qtthv_p

     ! Compute counter-gradient buoyancy flux
     thv2_p = max( 0., Beta_hl*hlthv_p + Beta_qt*qtthv_p )

     ! Restrict anisotropy by restricting buoyancy variance (NN09 Section 2.7)
!     if ( wden > 0.d0 ) then
     if ( Cw_high > 0.d0 ) then
        thv2_p = min( max( thv2_p, Cw_low ), Cw_high )
     else
        thv2_p = max( min( thv2_p, Cw_low ), Cw_high )
     end if
     wb_explicit = Lq*EH*goth002*thv2_p

     ! Update momentum stability function
     SM = SM25 + EM*L2*goth002*thv2_p
  else
     SM = SM25

     whl_explicit = 0.
     wqt_explicit = 0.
     wb_explicit  = 0.
  end if

  ! Compute diffusivities
  KM = Lq*SM
  KH = Lq*SH

  ! Compute counter-gradient fluxes of GEOS variables
  if ( MYNN_LEVEL >= 3 ) then
     wql_explicit = 0. ! temporary
     ws_explicit  = MAPL_CP*whl_explicit + MAPL_ALHL*wql_explicit
     wqv_explicit = wqt_explicit - wql_explicit
  end if

  !
  ! Compute total fluxes
  ! 

  if ( DOMF /= 0. .and. CONSISTENT_TYPE /= 0 ) then
     whl = -KH*dhldz + whl_explicit + whl_mf
     wqt = -KH*dqtdz + wqt_explicit + wqt_mf
  else
     whl = -KH*dhldz + whl_explicit
     wqt = -KH*dqtdz + wqt_explicit
  end if

  !
  ! Compute budget terms for second-order moments
  !

  ! Buoyancy produduction
  if ( domf /= 0. .and. consistent_type /= 0 ) then
     tket_B = -KH*N2 + wb_explicit + goth00*wthv_mf
  else
     tket_B = -KH*N2 + wb_explicit
  end if

  ! Mechanical production
  tket_M  = KM*S2
  hl2t_M  = -2.*whl*dhldz
  qt2t_M  = -2.*wqt*dqtdz
  hlqtt_M = -whl*dqtdz - wqt*dhldz

  ! Mass flux transport of TKE
  if ( CONSISTENT_TYPE == 0 ) then
     we    = -au*wu/( 1. - au )
     we_up = -au_up*wu_up/( 1. - au_up )
     
     tkee    = tke/( 1. - au )
     tkee_up = tke_up/( 1. - au_up )
     w2e_up  = twothirds*tkee_up
           
     tket_T_mf1 = ( Mu_up*tkee_up - Mu*tkee )/(dzle*rhoe)
     tket_T_mf2 = - max(0., E)*tkee/rhoe
     tket_T_mf3 =   max(0., D)*( we - wdet )**2./rhoe
     tket_T_mf4 = - ( 1. - au_up )*w2e_up*( we_up - we )/dzle

     tket_T_mf = tket_T_mf1 + tket_T_mf2 + tket_T_mf3 + tket_T_mf4
  else
     tket_T_mf1 = 0.
     tket_T_mf2 = 0.
     tket_T_mf3 = 0.
     tket_T_mf4 = 0.

     tket_T_mf = 0.
  end if

end subroutine mynn_tendency

!
!
!
subroutine implicit_M(IM, JM, LM, &                                                       ! in
                      th00, zlo, ple, u, v, h, qv, ql, Tv, tke, &                         ! in
                      Beta_hl, Beta_qt, L, qdiv, SM25, SH25, &                            ! in
                      ws_explicit, wqv_explicit, wql_explicit, whl_mf, wqt_mf, wthv_mf, & ! in
                      hl2, qt2, hlqt, &                                                   ! inout
                      tket_M, tket_B, hl2t_M, qt2t_M, hlqtt_M, &                          ! out
                      rhoe, dhldz, dqtdz, dqldz, S2, N2, &                                ! out
                      MYNN_LEVEL, DOMF, CONSISTENT_TYPE)                                  ! in

  integer, intent(in)                        :: IM, JM, LM, MYNN_LEVEL, CONSISTENT_TYPE
  real, intent(in)                           :: th00, DOMF
  real, dimension(IM,JM,LM), intent(in)      :: zlo, u, v, h, qv, ql, Tv 
  real, dimension(IM,JM,0:LM), intent(in)    :: ple, tke, Beta_hl, Beta_qt, L, qdiv, SM25, SH25, &
                                                ws_explicit, wqv_explicit, wql_explicit, whl_mf, wqt_mf, wthv_mf
  real, dimension(IM,JM,0:LM), intent(inout) :: hl2, qt2, hlqt
  real, dimension(IM,JM,0:LM), intent(out)   :: tket_M, tket_B, hl2t_M, qt2t_M, hlqtt_M, &
                                                dhldz, dqtdz, dqldz, S2, N2, rhoe

  integer :: i, j, k, kp1

  real                      :: idzlo, whl, wqt, whl_explicit, wqt_explicit, wb_explicit, goth00, &
                               T, Tve, KM, KH
  real, dimension(IM,JM,LM) :: hl, qt

  double precision :: Lq, L2 ! this precision may be unecessary

  goth00 = MAPL_GRAV/th00

  ! Compute conserved thermodynamic properties 
  do k = 1,LM
     do j = 1,JM
     do i = 1,IM
        hl(i,j,k) = h(i,j,k)/MAPL_CP - lvocp*ql(i,j,k)
        qt(i,j,k) = qv(i,j,k) + ql(i,j,k)
     end do
     end do
  end do

  do k = 1,LM-1

     kp1 = k + 1
     do j = 1,JM
     do i = 1,IM
        L2 = L(i,j,k)**2.
        Lq = real(L(i,j,k),8)*sqrt( real(2.*tke(i,j,k),8) )

        KM = Lq*SM25(i,j,k)
        KH = Lq*SH25(i,j,k)
        
        idzlo        = 1./( zlo(i,j,k) - zlo(i,j,kp1) )
        dhldz(i,j,k) = ( hl(i,j,k) - hl(i,j,kp1) )*idzlo
        dqtdz(i,j,k) = ( qt(i,j,k) - qt(i,j,kp1) )*idzlo
        dqldz(i,j,k) = ( ql(i,j,k) - ql(i,j,kp1) )*idzlo

        ! Compute density
        Tve         = 0.5*( Tv(i,j,k) + TV(i,j,kp1) )
        rhoe(i,j,k) = ple(i,j,k)/(MAPL_RGAS*Tve)

        N2(i,j,k) = goth00*( Beta_hl(i,j,k)*dhldz(i,j,k) + Beta_qt(i,j,k)*dqtdz(i,j,k) )
        S2(i,j,k) = (( u(i,j,k) - u(i,j,kp1) )*idzlo)**2. + ( (v(i,j,k) - v(i,j,kp1) )*idzlo)**2.

        whl_explicit = ws_explicit(i,j,k)/MAPL_CP - lvocp*wql_explicit(i,j,k)
        wqt_explicit = wqv_explicit(i,j,k) + wql_explicit(i,j,k)

        if ( mynn_level == 2 ) then
           wb_explicit = 0.

           ! Compute thermodyanamic (co-)variances from level-2.5 closure
           hl2(i,j,k)  = qdiv(i,j,k)*B2*L2*SH25(i,j,k)*dhldz(i,j,k)**2.
           qt2(i,j,k)  = qdiv(i,j,k)*B2*L2*SH25(i,j,k)*dqtdz(i,j,k)**2.
           hlqt(i,j,k) = qdiv(i,j,k)*B2*L2*SH25(i,j,k)*dhldz(i,j,k)*dqtdz(i,j,k)
        else
           wb_explicit = Beta_hl(i,j,k)*whl_explicit + Beta_qt(i,j,k)*wqt_explicit

           if ( mynn_level == 4 ) then
              ! Compute thermodyanamic (co-)variances from level-2.5 closure
              hl2(i,j,k)  = qdiv(i,j,k)*B2*L2*SH25(i,j,k)*dhldz(i,j,k)**2.
              hlqt(i,j,k) = qdiv(i,j,k)*B2*L2*SH25(i,j,k)*dhldz(i,j,k)*dqtdz(i,j,k)
           end if
        end if

        if ( DOMF /= 0. .and. CONSISTENT_TYPE /= 0 ) then
           tket_B(i,j,k) = -KH*N2(i,j,k)    + wb_explicit  + goth00*wthv_mf(i,j,k)
           whl           = -KH*dhldz(i,j,k) + whl_explicit + whl_mf(i,j,k)
           wqt           = -KH*dqtdz(i,j,k) + wqt_explicit + wqt_mf(i,j,k)
        else
           tket_B(i,j,k) = -KH*N2(i,j,k)    + wb_explicit
           whl           = -KH*dhldz(i,j,k) + whl_explicit
           wqt           = -KH*dqtdz(i,j,k) + wqt_explicit
        end if

        tket_M(i,j,k)  = KM*S2(i,j,k)
        hl2t_M(i,j,k)  = -2.*whl*dhldz(i,j,k)
        qt2t_M(i,j,k)  = -2.*wqt*dqtdz(i,j,k)
        hlqtt_M(i,j,k) = -whl*dqtdz(i,j,k) - wqt*dhldz(i,j,k)
     end do
     end do
  end do

end subroutine implicit_M

!
! mynn_length 
!
subroutine mynn_length(IM, JM, LM, &                                    ! in
                       local_flag, alpha1, alpha2, alpha3, alpha4, &    ! in
                       th00, ice_ramp, wb_surf, zle, zlo, q, N2, LMO, & ! in
                       thv, ple, pl, &                                  ! in
                       L, LS, LB, LT, w_star)                           ! out

  use MAPL_ConstantsMod, only: MAPL_KARMAN

  real, parameter :: qmin = 0.
  real, parameter :: zmax = 1.

  real, parameter :: cns = 2.7

  ! WRF-MYNN parameters
  real, parameter :: minzi = 300.
  real, parameter :: maxdz = 750.
  real, parameter :: mindz = 300.
  real, parameter :: zi2   = 10000.
  real, parameter :: h1    = min( maxdz, max( mindz, 0.3*zi2 ) )
!  real, parameter :: h2    = 0.5*h1

  real, parameter :: zi2ph1 = zi2 + h1


  integer, intent(in)                                  :: IM, JM, LM, local_flag
  real, intent(in)                                     :: th00, ice_ramp
  double precision, intent(in)                         :: alpha1, alpha2, alpha3, alpha4
  real, dimension(IM,JM), intent(in)                   :: wb_surf, LMO
  real, dimension(IM,JM,LM), intent(in)                :: zlo, thv, pl
  real, dimension(IM,JM,0:LM), intent(in)              :: zle, N2, ple
  double precision, dimension(IM,JM,0:LM), intent(in)  :: q
  double precision, dimension(IM,JM), intent(out)      :: w_star, LT
  double precision, dimension(IM,JM,0:LM), intent(out) :: L, LS, LB

  integer                            :: i, j, k, kk, kkm1
  real                               :: qdz, zeta
  double precision                   :: N, LF, LB_test, N_test
  double precision, dimension(IM,JM) :: work1, work2

  ! Compute integrals for turbulence length scale
  work1(:,:) = 1.d-5
  work2(:,:) = 1.d-5
  do k = 1,LM
     kk = LM - k + 1

     kkm1 = kk - 1
     do j = 1,JM
     do i = 1,IM
        if ( zle(i,j,kkm1) <= zi2ph1 ) then ! top limit of integration from WRF MYNN scheme
           qdz = max( 0.d0, ( 0.5*( q(i,j,kk) + q(i,j,kkm1) ) - qmin )*( zle(i,j,kkm1) - zle(i,j,kk) ) )

           work1(i,j) = work1(i,j) + zlo(i,j,kk)*qdz 
           work2(i,j) = work2(i,j) + qdz
        end if
     end do
     end do
  end do

  ! Compute turbulence length scale and CBL vertical velocity scale
  do j = 1,JM
  do i = 1,IM
     LT(i,j)     = alpha1*( work1(i,j)/work2(i,j) ) ! NN09 (54)
     w_star(i,j) = ( max( 0.d0, wb_surf(i,j) )*LT(i,j) )**onethird
  end do
  end do

  ! Compute master length scale
  do k = 1,LM-1
     do j = 1,JM
     do i = 1,IM
        ! Compute non-dimensional height
        zeta = zle(i,j,k)/LMO(i,j)

        ! WRF MYNN version
        if ( zeta > 0. ) then
           LS(i,j,k) = MAPL_KARMAN*zle(i,j,k)/( 1. + cns*min(zmax, zeta) )
        else
           LS(i,j,k) = MAPL_KARMAN*zle(i,j,k)*( 1. - alpha4*zeta )**0.2
        end if
        
        ! Compute buoyancy length scale
        if ( N2(i,j,k) > 0. ) then
           N = sqrt( N2(i,j,k) )

           if ( local_flag == 0 ) then
              call boulac(IM, JM, LM, i, j, k, &                    ! in
                          th00, ice_ramp, zlo, zle, pl, ple, thv, & ! in
                          q(i,j,k), N, &                            ! in
                          LB_test)                                  ! out    

              N_test    = q(i,j,k)/LB_test
              LB(i,j,k) = alpha2*( LB_test )*( 1.d0 + alpha3/alpha2*sqrt( w_star(i,j)/( N_test*LT(i,j) ) ) )
              LF        = alpha2*( LB_test )
           else
              LB(i,j,k) = alpha2*( q(i,j,k)/N )*( 1.d0 + alpha3/alpha2*sqrt( w_star(i,j)/( N*LT(i,j) ) ) )
              LF        = alpha2*( q(i,j,k)/N )              
           end if
        else
           LB(i,j,k) = 1.d+10
           LF        = LB(i,j,k)
        end if
        
        ! Harmonically average length scales
        L(i,j,k) = min( LF, LB(i,j,k)/( LB(i,j,k)/LS(i,j,k) + LB(i,j,k)/LT(i,j) + 1.d0 ) ) ! NN09 (52)
     end do
     end do
  end do

  LS(:,:,LM) = 0.
  L(:,:,LM)  = L(:,:,LM-1)

end subroutine mynn_length

!
! initialize_mynn
!
subroutine initialize_mynn(IM, JM, LM, &
                           local_flag, alpha1, alpha2, alpha3, alpha4, &
                           th00, ice_ramp, hl, qt, tke, hl2, qt2, hlqt, q, &
                           zle, zlo, S2, N2, &
                           u_star, wb_surf, LMO, &
                           thv, ple, pl)

use MAPL_ConstantsMod, only: MAPL_KARMAN

integer, parameter :: niter = 5

real, parameter :: pmz = 1.
real, parameter :: phh = 1.
real, parameter :: flt = 0.
real, parameter :: flq = 0.

real, parameter :: phm = phh*B2/(B1*pmz)**twothirds
       
integer, intent(in)                                    :: IM, JM, LM, local_flag
real, intent(in)                                       :: th00, ice_ramp
double precision, intent(in)                           :: alpha1, alpha2, alpha3, alpha4
real, dimension(IM,JM), intent(in)                     :: u_star, wb_surf, LMO
real, dimension(IM,JM,0:LM), intent(in)                :: zle, S2, N2, ple
real, dimension(IM,JM,LM), intent(in)                  :: zlo, hl, qt, thv, pl
real, dimension(IM,JM,0:LM), intent(inout)             :: tke, hl2, qt2, hlqt
double precision, dimension(IM,JM,0:LM), intent(inout) :: q

integer                                 :: iter, i, j, k, kp1
real                                    :: idzlo, L2
double precision, dimension(IM,JM)      :: w_star
real, dimension(IM,JM,0:LM)             :: SM2, SH2, dhldz, dqtdz
double precision, dimension(IM,JM)      :: LT
double precision, dimension(IM,JM,0:LM) :: GM, GH, L, LS, LB


!
do k = 1,LM-1
   kp1 = k + 1

   do j = 1,JM
   do i = 1,IM
      idzlo        = 1./ ( zlo(i,j,k) - zlo(i,j,kp1) )
      dhldz(i,j,k) = ( hl(i,j,k) - hl(i,j,kp1) )*idzlo
      dqtdz(i,j,k) = ( qt(i,j,k) - qt(i,j,kp1) )*idzlo

      GH(i,j,k) = -N2(i,j,k) ! This is actually GH divided by L2/q2
      GM(i,j,k) = S2(i,j,k)  ! "    "  "        GM "       "  "

      call mynn_l2(GM(i,j,k), GH(i,j,k), SM2(i,j,k), SH2(i,j,k))
   end do
   end do
end do

! First pass of TKE initialization
do j = 1,JM
do i = 1,IM
   tke(i,j,LM-1) = 0.5*u_star(i,j)**2.*(B1*pmz)**twothirds
   q(i,j,LM-1)   = sqrt( max( 1.d-10, real(2.*tke(i,j,LM-1),8) ) )
end do
end do

do k = 1,LM-2
   do j = 1,JM
   do i = 1,IM
      tke(i,j,k) = 0.
      q(i,j,k)   = 1.d-10
   end do
   end do
end do

! Iterate to initialize TKE
do iter = 1,niter
   call mynn_length(IM, JM, LM, &                                    ! in
                    local_flag, alpha1, alpha2, alpha3, alpha4, &    ! in
                    th00, ice_ramp, wb_surf, zle, zlo, q, N2, LMO, & ! in
                    thv, ple, pl, &                                  ! in
                    L, LS, LB, LT, w_star)                           ! out      

   do k = 1,LM-1
      do j = 1,JM
      do i = 1,IM
         L2 = L(i,j,k)**2.d0

         tke(i,j,k) = 0.5*B1*L2*( SM2(i,j,k)*GM(i,j,k) + SH2(i,j,k)*GH(i,j,k) )
         q(i,j,k)   = sqrt( max( 1.d-10, real(2.*tke(i,j,k),8) ) )
      end do
      end do
   end do
end do

! Initialize second-order thermodynamic moments
do k = 1,LM-1

   do j = 1,JM
   do i = 1,IM
      L2 = L(i,j,k)**2.d0

      hl2(i,j,k)  = B2*L2*SH2(i,j,k)*dhldz(i,j,k)**2.
      qt2(i,j,k)  = B2*L2*SH2(i,j,k)*dqtdz(i,j,k)**2.
      hlqt(i,j,k) = B2*L2*SH2(i,j,k)*dhldz(i,j,k)*dqtdz(i,j,k)
   end do
   end do
end do

end subroutine initialize_mynn

!
! mynn_l2
!
subroutine mynn_l2(GM, GH, SM2, SH2)

  double precision, intent(in) :: GM, GH
  real, intent(out)            :: SM2, SH2

  double precision :: Ri, Rf

  Ri  = -GH/max( GM, 1.d-10 )
  Rf  = min( Rfc, Ri1*( Ri + Ri2 - sqrt( Ri**2.d0 - Ri3*Ri + Ri4 ) ) ) ! NN09 (A11)

  SH2 = SHc*( ( Rfc - Rf )/( 1.d0 - Rf ) )     ! NN09 (A4)
  SM2 = SMc*( ( Rf1 - Rf )/( Rf2 - Rf )  )*SH2 ! NN09 (A3)

end subroutine mynn_l2

!
! mynn_l25
!
subroutine mynn_l25(L, q, tke, dhldz, dqtdz, S2, N2, &    ! in
                    SM25, SH25, &                         ! out (MYNN level 2.5 state)
                    Phi1, Phi2, Phi3, Phi4, Phi5, D_25, & ! out (needed for MYNN level-3)
                    L2, GH, L2GH, q2, qdiv, qdiv2)        ! out (needed for MYNN level-3)

  double precision, intent(in)  :: L, q
  real, intent(in)              :: tke, dhldz, dqtdz, S2, N2
  real, intent(out)             :: SM25, SH25
  double precision, intent(out) :: Phi1, Phi2, Phi3, Phi4, Phi5, D_25, &
                                   L2, GH, L2GH, q2, qdiv, qdiv2

  real :: SM2, SH2

  double precision :: GM, L2GM, q22

  ! Compute some intermediate quantities
  L2   = L**2.d0
  q2   = q**2.d0
  GH   = -N2 ! This is actually GH divided by L2/q2
  GM   = S2  ! "    "  "        GM "       "  "
  L2GH = L2*GH
  L2GM = L2*GM

  ! Compute Level-2 closure quantities
  call mynn_l2(GM, GH, SM2, SH2)
  q22 = B1*L2*( SM2*GM + SH2*GH ) ! NN09 (A1))

  ! Compute (1-alpha) and (1-alpha)^2 factors from HL88
  if ( q2 < q22 ) then
     qdiv2 = q2/q22
     qdiv  = sqrt(qdiv2)
  else
     qdiv2 = 1.d0
     qdiv  = 1.d0
  end if

  ! Compute useful intermediate quantities
  Phi1 = q2   - e1c*L2GH*qdiv2                ! NN09 (33)
  Phi2 = q2   - e2c*L2GH*qdiv2                ! NN09 (34)
  Phi3 = Phi1 + e3c*L2GH*qdiv2                ! NN09 (35)
  Phi4 = Phi1 - e4c*L2GH*qdiv2                ! NN09 (36)
  Phi5 =        e5c*L2GM*qdiv2                ! NN09 (37)
  D_25 = max( 1.d-20, Phi2*Phi4 + Phi3*Phi5 ) ! NN09 (31)

  ! Compute stability functions
  if ( q2 < q22 ) then
     SM25 = qdiv*SM2
     SH25 = qdiv*SH2
  else
     SM25 = q2*A1*( ( Phi3 - 3.d0*C1*Phi4 )/D_25 ) ! NN09 (27)
     SH25 = q2*A2*( ( Phi2 + 3.d0*C1*Phi5 )/D_25 ) ! NN09 (28)
  end if

end subroutine mynn_l25

!
! mynn_l3
!
subroutine mynn_l3(th00, &                                                    ! in 
                   Phi1, Phi3, Phi4, Phi5, D_25, L2, GH, L2GH, qdiv, qdiv2, & ! in (MYNN level-2.5 state)
                   Phi2, q2, &                                                ! inout         
                   EM, EH, Cw_low, Cw_high)                                   ! out (MYNN level-3 state)

  real, intent(in)                :: th00
  double precision, intent(in)    :: Phi1, Phi3, Phi4, Phi5, D_25, L2, GH, L2GH, qdiv, qdiv2
  double precision, intent(inout) :: Phi2, q2
  real, intent(out)               :: EM, EH, Cw_low, Cw_high

  real             :: Cw_25, goth002
  double precision :: D_p, wden

  goth002 = (MAPL_GRAV/th00)**2.

  ! Limit q2 so that L/q is less than 1/N for N2 > 0 (NN09 Section 2.7)
  if ( q2/L2 < -GH ) then
     q2 = -L2GH ! NN09 (56)
  end if

  !
  Phi2 = q2 - e2c*L2GH*qdiv2                                                  ! NN09 (34)
  D_p  = max( 1.d-20, Phi2*( Phi4 - Phi1 + q2 ) + Phi5*( Phi3 - Phi1 + q2 ) ) ! NN09 (32)

  ! Compute level-3 momentum stability function
  EH = qdiv*eHc*( ( Phi2 + Phi5 )/D_p )        ! NN09 (48)
  EM = qdiv*eMc*( ( Phi3 - Phi4 )/(D_p*L2GH) ) ! NN09 (47)

  !
  wden = ( 1.d0 - C3 )*goth002*L2*qdiv2*( e4c*Phi2 - e3c*Phi5 )
  if ( wden /= 0.d0 ) then
     Cw_25   = Phi1*( Phi2 + 3.d0*C1*Phi5 )/(3.d0*D_25) ! NN09 (59)
     Cw_low  = q2*( 0.12 - Cw_25 )*( D_p/wden )
     Cw_high = q2*( 0.76 - Cw_25 )*( D_p/wden )
  end if

end subroutine mynn_l3

!
! boulac
!
subroutine boulac(IM, JM, LM, i, j, kz, &                  ! in
                  th00, ice_ramp, zl, zle, pl, ple, thv, & ! in
                  q, N, &                                  ! in
                  LB)                                      ! out

  use edmf_mod, only: condensation_edmfA, condensation_edmf

  implicit none

  integer, intent(in)                     :: IM, JM, LM, i, j, kz
  real, intent(in)                        :: th00, ice_ramp
  double precision, intent(in)            :: q, N
  real, dimension(IM,JM,LM), intent(in)   :: zl, pl, thv
  real, dimension(IM,JM,0:LM), intent(in) :: zle, ple
  double precision, intent(out)           :: LB

  integer :: k, km1, kp1
  real    :: thvp, lup, ldown, w2p, w2p_next, ifac, &
             lup0, ldown0, dz, wf, goth00, B

  goth00 = mapl_grav/th00

  lup0   = zl(i,j,kz)  - zle(i,j,kz)
  ldown0 = zle(i,j,kz) - zl(i,j,kz+1)

  ! If LB smaller that half of dz, then use local bouyancy length scale
  if ( q/N <= min(lup0,ldown0) ) then
     LB = q/N  
     return
  end if

  ! Initialize parcel
  ifac = ( zle(i,j,kz) - zl(i,j,kz) )/( zl(i,j,kz) - zl(i,j,kz+1) )
  thvp = thv(i,j,kz+1) + ifac*( thv(i,j,kz) - thv(i,j,kz+1) )

  ! up
  k   = kz
  lup = lup0
  w2p = q**2 - (lup0*N)**2.
up:  do while ( w2p > 0. .and. k >= 2 )
     km1 = k - 1

     dz = zl(i,j,km1) - zl(i,j,k)

     B = goth00*( thvp - 0.5*( thv(i,j,k) + thv(i,j,km1) ) )

     w2p_next = w2p + 2.*dz*B

     if ( w2p_next > 0. ) then
        lup = lup + dz
        w2p = w2p_next
        k   = k - 1
     else
        lup = lup - 0.5*w2p/B

        exit up
     end if
  end do up

  ! down
  k     = kz + 1
  ldown = ldown0
  w2p   = -q**2. + (ldown0*N)**2.
down:  do while ( w2p < 0. .and. k <= LM-1 )
     kp1 = k + 1

     dz = zl(i,j,k) - zl(i,j,kp1)

     B = goth00*( thvp - 0.5*( thv(i,j,k) + thv(i,j,kp1) ) )

     w2p_next = w2p + 2.*dz*B

     if ( w2p_next < 0. ) then
        ldown = ldown + dz

        ! Check proximity to surface
        if ( k <= LM-2 ) then 
           w2p = w2p_next
           k   = k + 1
        else
           ldown = ldown + zl(i,j,LM)

           exit down
        end if
     else
        ldown = ldown - 0.5*w2p/B

        exit down
     end if
  end do down

  ! Geometrically average each length scale
  LB = min( q/N, min( zle(i,j,kz), sqrt(lup*ldown) ) )

end subroutine boulac

!
! mynn_predict_correct
!
subroutine mynn_predict_correct(IM, JM, LM, &                                         ! in
                                mynn_level, wrf_cg_flag, domf, consistent_type, &     ! in
                                th00, zle, ple, rhoe, tke, hl2, qt2, hlqt, &          ! in
                                dhldz, dqtdz, dqldz, N2, S2, &                        ! in (gradient information)
                                Beta_hl, Beta_qt, qdiv, L, SM25, SH25, EM, EH, &      ! in
                                au, Mu, E, D, wu, wdet, &                             ! in (for consistent partitioning of TKE)
                                whl_mf, wqt_mf, wthv_mf, &                            ! in (for non-consistent partitioning of TKE)
                                K_tke, itau, &                                        ! inout
                                tket_M, tket_B, tket_T_mf, hl2t_M, qt2t_M, hlqtt_M, & ! inout
                                tket_T_mf1, tket_T_mf2, tket_T_mf3, tket_T_mf4)       ! inout

  implicit none

  integer, intent(in)                     :: IM, JM, LM, mynn_level, consistent_type, wrf_cg_flag
  real, intent(in)                        :: domf, th00
  real, dimension(IM,JM,0:LM), intent(in) :: zle, ple, rhoe, Beta_hl, Beta_qt, &
                                             dhldz, dqtdz, dqldz, N2, S2, &
                                             L, qdiv, SM25, SH25, EM, EH, &
                                             whl_mf, wqt_mf, wthv_mf, tke, hl2, qt2, hlqt, &
                                             au, Mu, E, D, wu, wdet
  real, dimension(IM,JM,LM), intent(out)     :: K_tke
  real, dimension(IM,JM,0:LM), intent(inout) :: itau, &
                                                tket_M, tket_B, tket_T_mf, hl2t_M, qt2t_M, hlqtt_M, &
                                                tket_T_mf1, tket_T_mf2, tket_T_mf3, tket_T_mf4

  integer :: i, j, k, km1, kp1
  real :: KH, hl2_25, qt2_25, hlqt_25, ws_explicit, wqv_explicit, wql_explicit, &
          tket_M_test, tket_B_test, tket_T_mf_test, hl2t_M_test, qt2t_M_test, hlqtt_M_test, &
          tket_T_mf1_test, tket_T_mf2_test, tket_T_mf3_test, tket_T_mf4_test, dzle
  double precision :: L2, Lq, q, L_double
  real, dimension(IM,JM,0:LM) :: KM

  do k = 1,LM-1

     km1 = k - 1
     kp1 = k + 1
     do j = 1,JM
     do i = 1,IM
        dzle = zle(i,j,km1) - zle(i,j,k)

        q        = sqrt( max( 1.d-10, real(2.*tke(i,j,k),8) ) )
        L_double = real( L(i,j,k), 8)

        itau(i,j,k) = q/L_double
        L2          = L_double**2.d0
        Lq          = L_double*real( sqrt(2.*tke(i,j,k)), 8 )

        ! Compute thermodyanamic (co-)variances from level-2.5 closure
        hl2_25  = qdiv(i,j,k)*B2*L2*SH25(i,j,k)*dhldz(i,j,k)**2.
        qt2_25  = qdiv(i,j,k)*B2*L2*SH25(i,j,k)*dqtdz(i,j,k)**2.
        hlqt_25 = qdiv(i,j,k)*B2*L2*SH25(i,j,k)*dhldz(i,j,k)*dqtdz(i,j,k)

        ! Compute diffusivities and second-order tendencies
        call mynn_tendency(mynn_level, wrf_cg_flag, domf, consistent_type, &                                               ! in
                           th00, rhoe(i,j,k), Beta_hl(i,j,k), Beta_qt(i,j,k), hl2(i,j,k), qt2(i,j,k), hlqt(i,j,k), &       ! in
                           dzle, dhldz(i,j,k), dqtdz(i,j,k), dqldz(i,j,k), N2(i,j,k), S2(i,j,k), &                         ! in
                           L2, Lq, SM25(i,j,k), SH25(i,j,k), EM(i,j,k), EH(i,j,k), hl2_25, qt2_25, hlqt_25, &              ! in
                           tke(i,j,k), tke(i,j,km1), Mu(i,j,k), Mu(i,j,km1), au(i,j,k), au(i,j,km1), &                     ! in (for tket_T_mf)
                           wu(i,j,k), wu(i,j,km1), E(i,j,k), D(i,j,k), wdet(i,j,k), &                                      ! in (for tket_T_mf)
                           whl_mf(i,j,k), wqt_mf(i,j,k), wthv_mf(i,j,k), &                                                 ! in
                           KM(i,j,k), KH, ws_explicit, wqv_explicit, wql_explicit, &                                 ! out
                           tket_M_test, tket_B_test, tket_T_mf_test, hl2t_M_test, qt2t_M_test, hlqtt_M_test, &       ! out
                           tket_T_mf1_test, tket_T_mf2_test, tket_T_mf3_test, tket_T_mf4_test)                       ! out

        ! Improved Euler's method
        tket_M(i,j,k)     = 0.5*( tket_M(i,j,k)     + tket_M_test )
        tket_B(i,j,k)     = 0.5*( tket_B(i,j,k)     + tket_B_test )
        tket_T_mf(i,j,k)  = 0.5*( tket_T_mf(i,j,k)  + tket_T_mf_test )
        hl2t_M(i,j,k)     = 0.5*( hl2t_M(i,j,k)     + hl2t_M_test )
        qt2t_M(i,j,k)     = 0.5*( qt2t_M(i,j,k)     + qt2t_M_test )
        hlqtt_M(i,j,k)    = 0.5*( hlqtt_M(i,j,k)    + hlqtt_M_test )

        tket_T_mf1(i,j,k) = 0.5*( tket_T_mf1(i,j,k) + tket_T_mf1_test )
        tket_T_mf2(i,j,k) = 0.5*( tket_T_mf2(i,j,k) + tket_T_mf2_test )
        tket_T_mf3(i,j,k) = 0.5*( tket_T_mf3(i,j,k) + tket_T_mf3_test )
        tket_T_mf4(i,j,k) = 0.5*( tket_T_mf4(i,j,k) + tket_T_mf4_test )
     end do
     end do
  end do

  ! Compute TKE diffusivity
  do k = 2,LM-1
     km1 = k - 1
     do j = 1,JM
     do i = 1,IM
        K_tke(i,j,k) = 1.5*( KM(i,j,k) + KM(i,j,km1) )
     end do
     end do
  end do
  K_tke(:,:,1)  = K_tke(:,:,2)
  K_tke(:,:,LM) = K_tke(:,:,LM-1)

end subroutine mynn_predict_correct



!
! entrain_mynn
!
subroutine entrain_mynn(IM, JM, LM, &
                        th00, zl, zle, omega, thl, qt, thv, ac, &
                        thlv, ktop, zi, gamma_ml, gamma_fa)

integer, intent(in)                     :: IM, JM, LM
real, intent(in)                        :: th00
real, dimension(IM,JM,LM), intent(in)   :: zl, thl, qt, thv, ac
real, dimension(IM,JM,0:LM), intent(in) :: zle, omega
real, dimension(IM,JM,LM), intent(out)  :: thlv
integer, intent(out), dimension(IM,JM)  :: ktop
real, intent(out),dimension(IM,JM)      :: zi, gamma_ml, gamma_fa    

integer                   :: i, j, k, kflip
real                      :: a, b, c, thlv_ml, thlv_fa, we, w_ml
real                      :: goth00
logical, dimension(IM,JM) :: search_flag

goth00 = MAPL_GRAV/th00

! Find cloud top
ktop(:,:) = -1
search_flag = .true.
do k = 1,LM
   kflip = LM - k + 1
   do j = 1,JM
   do i = 1,IM
      if ( search_flag(i,j) .and. ac(i,j,kflip) > 0.01 ) then
         ktop(i,j)        = kflip
         search_flag(i,j) = .false.
      end if
   end do
   end do
end do

! Entrainment closure
thlv = thl*( 1. + 0.622*qt )
do j = 1,JM
do i = 1,IM
   gamma_ml(i,j) = ( thlv(i,j,ktop(i,j)+1) - thlv(i,j,ktop(i,j)+2) )/( zl(i,j,ktop(i,j)+1) - zl(i,j,ktop(i,j)+2) )
   gamma_fa(i,j) = ( thlv(i,j,ktop(i,j)-2) - thlv(i,j,ktop(i,j)-1) )/( zl(i,j,ktop(i,j)-2) - zl(i,j,ktop(i,j)-1) )

   a = 0.5*( gamma_fa(i,j) - gamma_ml(i,j) )
   b = - ( thlv(i,j,ktop(i,j)-1) - gamma_fa(i,j)*( zl(i,j,ktop(i,j)-1) -  zle(i,j,ktop(i,j)) ) ) &
        + ( thlv(i,j,ktop(i,j)+1) + gamma_ml(i,j)*(  zle(i,j,ktop(i,j)) - zl(i,j,ktop(i,j)+1) ) )
   c =  ( zle(i,j,ktop(i,j)) - zle(i,j,ktop(i,j)+1) )&
        *( thlv(i,j,ktop(i,j)) - ( thlv(i,j,ktop(i,j)+1) + gamma_ml(i,j)*( zl(i,j,ktop(i,j)) - zl(i,j,ktop(i,j)+1) ) ) )

   if ( b > 0. ) then
      b = thlv(i,j,ktop(i,j)+1) - thlv(i,j,ktop(i,j)-1)
      c = ( zle(i,j,ktop(i,j)) - zle(i,j,ktop(i,j)+1) )*( thlv(i,j,ktop(i,j)) - thlv(i,j,ktop(i,j)+1) )

      zi(i,j) = zle(i,j,ktop(i,j)) - c/b
   else
      zi(i,j) = zle(i,j,ktop(i,j)) - ( -b - sqrt( b**2. - 4.*a*c  ) )/( 2.*a )
   end if

   thlv_ml = thlv(i,j,ktop(i,j)+1) + gamma_ml(i,j)*(             zi(i,j) - zl(i,j,ktop(i,j)+1) )
   thlv_fa = thlv(i,j,ktop(i,j)-1) - gamma_fa(i,j)*( zl(i,j,ktop(i,j)-1) -             zi(i,j) )

!   we = 0.16*tke(i,j,ktop(i,j))**1.5/(L(i,j,ktop(i,j))*goth00*( thlv_ml - thlv_fa ))
   
   w_ml = (MAPL_RGAS/MAPL_GRAV)*omega(i,j,k)*( thv(i,j,ktop(i,j)) + thv(i,j,ktop(i,j)) )

end do
end do

end subroutine entrain_mynn

end module mynn

!!$     if ( WQL_TYPE == 0 ) then
!!$        wql_explicit = 0.
!!$     else
!!$        ifac    = (zle - zlo(i,j,k+1))*idzlo
!!$        ac_half = ac(i,j,k+1) + ifac*(ac - ac(i,j,k+1))
!!$        wql     = ac_half*(  A*( -KH*dqtdz(i,j,k) + wqt_explicit ) &
!!$                           - B*( -KH*dhldz(i,j,k) + whl_explicit ) )
!!$              
!!$        wql_explicit = wql + KH*dqldz
!!$     end if
