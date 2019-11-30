module mynn

use MAPL_ConstantsMod, only: MAPL_GRAV, MAPL_CP, MAPL_ALHL, MAPL_RGAS, MAPL_RVAP

implicit none

double precision, parameter :: onethird = 1./3.

double precision, parameter :: th00 = 300.

double precision, parameter :: gocp    = MAPL_GRAV/MAPL_CP
double precision, parameter :: lvocp   = MAPL_ALHL/MAPL_CP
double precision, parameter :: kappa   = MAPL_RGAS/MAPL_CP
double precision, parameter :: goth00  = MAPL_GRAV/th00
double precision, parameter :: ep2     = MAPL_RVAP/MAPL_RGAS - 1.
double precision, parameter :: goth002 = goth00**2.

! MYNN constants 
double precision, parameter :: Pr     = 0.74
double precision, parameter :: gamma1 = 0.235
double precision, parameter :: B1     = 24.
double precision, parameter :: B2     = 15.
double precision, parameter :: C2     = 0.729
double precision, parameter :: C3     = 0.34
double precision, parameter :: C4     = 0.
double precision, parameter :: C5     = 0.2

double precision, parameter :: A1     = B1*( 1. - 3.*gamma1 )/6.
double precision, parameter :: C1     = gamma1 - 1./(3.*A1*B1**onethird)
double precision, parameter :: A2     = A1*( gamma1 - C1 )/(gamma1*Pr)
double precision, parameter :: gamma2 = B2/B1*( 1. - C3 ) + 2.*A1/B1*( 3. - 2.*C2 )

double precision, parameter :: e1c = 3.*A2*B2*( 1. - C3 )
double precision, parameter :: e2c = 9.*A1*A2*( 1. - C2 )
double precision, parameter :: e3c = 9.*A2**2.*( 1. - C2 )*( 1. - C5 )
double precision, parameter :: e4c = 12.*A1*A2*( 1. - C2 )
double precision, parameter :: e5c = 6.*A1**2.
double precision, parameter :: eMc = 3.*A1*( 1. - C3 )
double precision, parameter :: eHc = 3.*A2*( 1. - C3 )

double precision, parameter :: Rfc = gamma1/( gamma1 + gamma2 )                                                ! NN09 (A10)
double precision, parameter :: F1  = B1*( gamma1 - C1 ) + 2.*A1*( 3. - 2.*C2 ) + 3.*A2*( 1. - C2 )*( 1. - C5 ) ! NN09 (A6)
double precision, parameter :: F2  = B1*( gamma1 + gamma2 ) - 3.*A1*( 1. - C2 )                                ! NN09 (A7)
double precision, parameter :: Rf1 = B1*( gamma1 - C1 )/F1                                                     ! NN09 (A8)
double precision, parameter :: Rf2 = B1*gamma1/F2                                                              ! NN09 (A9)

double precision, parameter :: Smc = A1/A2*F1/F2 
double precision, parameter :: Shc = 3.*A2*( gamma1 + gamma2 )

double precision, parameter :: Ri1 = 0.5/Smc             ! NN09 (A12)
double precision, parameter :: Ri2 = Rf1*Smc             ! NN09 (A13)
double precision, parameter :: Ri3 = 4.*Rf2*Smc - 2.*Ri2 ! NN09 (A14)
double precision, parameter :: Ri4 = Ri2**2.

! Test
logical :: initialized_mynn = .false.

contains

!
!
!
subroutine run_mynn(IM, JM, LM, &
                    ple, zle, zlo, &
                    u, v, T, qv, ql, qi, tke, hl2, qt2, hlqt, &
                    u_star, H, E, thv, ac, whl_mf, wqt_mf, wthv_mf, &
                    Kh, Km, K_tke, K_tpe, itau, ws_cg, wqv_cg, wql_cg, &
                    tket_M, tket_B, hl2t_M, qt2t_M, hlqtt_M, &
                    DOMF, MYNN_LEVEL)

  use MAPL_ConstantsMod, only: MAPL_KARMAN
  use MAPL_SatVaporMod, only: MAPL_EQsat

  integer, intent(in)                        :: IM, JM, LM, MYNN_LEVEL
  real, intent(in)                           :: DOMF
  real, dimension(IM,JM), intent(in)         :: u_star, H, E
  real, dimension(IM,JM,LM), intent(in)      :: zlo, u, v, T, qv, ql, qi, ac, thv
  real, dimension(IM,JM,0:LM), intent(in)    :: ple, zle, hl2, qt2, hlqt, whl_mf, wqt_mf, wthv_mf
  real, dimension(IM,JM,0:LM), intent(inout) :: tke
  real, dimension(IM,JM,0:LM), intent(out)   :: Kh, Km, itau, ws_cg, wqv_cg, wql_cg, &
                                                tket_M, tket_B, hl2t_M, qt2t_M, hlqtt_M
  real, dimension(IM,JM,LM), intent(out)     :: K_tke, K_tpe

  integer :: i, j, k, kp1, km1
  double precision :: GH, GM, dhldz, dqtdz, dqldz, idzlo, exner, ifac, iexner, &
                      Sm2, Sh2, Sm, Sh, Cw_low, Cw_high, wrk1, &
                      Cw_25, rho_surf, LMO, whl, wqt, Ri, Rf, &
                      whl_cg, wqt_cg, wb_cg, Lq, wql, &
                      ac_half, T_half, ql_half, Tl, qs, dqs, &
                      q2, q22, EM, EH, Phi1, Phi2, Phi3, Phi4, Phi5, &
                      D_25, D_p, wden, qdiv, qdiv2, L2, L2GM, L2GH, &
                      hl2_25, qt2_25, hlqt_25, hlthv, qtthv, thv2, &
                      hlthv_25, qtthv_25, thv2_25, hlthv_p, qtthv_p, thv2_p

  double precision, dimension(IM,JM)      :: wb_surf
  double precision, dimension(IM,JM,LM)   :: hl, qt
  double precision, dimension(IM,JM,0:LM) :: S2, N2, zeta, A, B, Beta_hl, Beta_qt, L, q

  ! Test
  double precision :: w2, wb, tau, whl_test, wqt_test, wb_test, Kh_test

  ! Test
  if (.not. initialized_mynn) then
     tke(:,:,LM-1) = 0.1
     initialized_mynn = .true.
  end if

  ! Compute conserved thermodynamic properties
  do k = 1,LM
     do j = 1,JM
     do i = 1,IM
        hl(i,j,k) = T(i,j,k) + gocp*zlo(i,j,k) - lvocp*ql(i,j,k)
        qt(i,j,k) = qv(i,j,k) + ql(i,j,k)
     end do
     end do
  end do

  ! Compute shear and buoyancy frequencies
  do k = 1,LM-1

     kp1 = k + 1
     do j = 1,JM
     do i = 1,IM
        q(i,j,k) = sqrt(max(1.E-10, 2.*tke(i,j,k)))

        idzlo = 1./( zlo(i,j,k) - zlo(i,j,kp1) )
        dhldz = ( hl(i,j,k) - hl(i,j,kp1) )*idzlo
        dqtdz = ( qt(i,j,k) - qt(i,j,kp1) )*idzlo
        dqldz = ( ql(i,j,k) - ql(i,j,kp1) )*idzlo

        exner  = (ple(i,j,k)*1.E-5)**kappa
        ifac   = (zle(i,j,k) - zlo(i,j,kp1))*idzlo
        iexner = 1./exner

        ac_half = ac(i,j,kp1) + ifac*(ac(i,j,k) - ac(i,j,kp1))
        T_half  = T(i,j,kp1)  + ifac*(T(i,j,k)  - T(i,j,kp1))
        ql_half = ql(i,j,kp1) + ifac*(ql(i,j,k) - ql(i,j,kp1))

        Tl = T_half - lvocp*ql_half
        qs = MAPL_EQsat(Tl, real(ple(i,j,k), 8), dqs)

        wrk1 = lvocp*iexner - (1. + ep2)*th00
        A(i,j,k) = 1./( 1. + lvocp*dqs )
        B(i,j,k) = A(i,j,k)*exner*dqs

        Beta_hl(i,j,k) = iexner*( 1. - wrk1*ac_half*B(i,j,k) )
        Beta_qt(i,j,k)  = ep2*th00 + wrk1*ac_half*A(i,j,k)

        N2(i,j,k) = goth00*( Beta_hl(i,j,k)*dhldz + Beta_qt(i,j,k)*dqtdz )
        S2(i,j,k) = (( u(i,j,k) - u(i,j,kp1) )*idzlo)**2. + ( (v(i,j,k) - v(i,j,kp1) )*idzlo)**2.
     end do
     end do
  end do

  !
  do j = 1,JM
  do i = 1,IM
     rho_surf     = ple(i,j,LM)/(MAPL_RGAS*thv(i,j,LM))
     wb_surf(i,j) = goth00*( H(i,j)/(MAPL_CP*rho_surf) + ep2*th00*E(i,j)/rho_surf )
     LMO          = -u_star(i,j)**3./(MAPL_KARMAN*wb_surf(i,j))

     tket_B(i,j,LM) = wb_surf(i,j)
  end do
  end do

  call mynn_length(IM, JM, LM, wb_surf, zle, zlo, q, N2, LMO, L)

  do k = 1,LM-1

     kp1 = k + 1
     do j = 1,JM
     do i = 1,IM
        !
        idzlo   = 1./( zlo(i,j,k) - zlo(i,j,kp1) )
        dhldz   = (hl(i,j,k) - hl(i,j,kp1))*idzlo
        dqtdz   = (qt(i,j,k) - qt(i,j,kp1))*idzlo
        dqldz   = (ql(i,j,k) - ql(i,j,kp1))*idzlo

        !
        L2   = L(i,j,k)**2.
        GH   = -N2(i,j,k) ! This is actually GH divided by L2/q2
        GM   = S2(i,j,k)  ! "    "  "        GM "       "  "
        L2GH = L2*GH
        L2GM = L2*GM
        q2   = q(i,j,k)**2.

        ! Compute Level-2 closure quantities
        Ri  = -GH/max( GM, 1.E-10 )
        Rf  = min( Rfc, Ri1*( Ri + Ri2 - sqrt( Ri**2. - Ri3*Ri + Ri4 ) ) ) ! NN09 (A11)
        Sh2 = Shc*( Rfc - Rf )/( 1. - Rf )                                 ! NN09 (A4)
        Sm2 = Smc*( Rf1 - Rf )/( Rf2 - Rf )*Sh                             ! NN09 (A3)
        q22 = B1*L2*( Sm2*GM + Sh2*GH )                                    ! NN09 (A1))

        ! Compute (1-alpha) and (1-alpha)^2 factors from HL88
        if (q2 < q22) then
           qdiv2 = q2/q22
           qdiv  = sqrt(qdiv2)
        else
           qdiv2 = 1.
           qdiv  = 1.
        end if

        ! Compute useful intermediate quantities
        Phi1 = q2   - e1c*L2GH*qdiv2                ! NN09 (33)
        Phi2 = q2   - e2c*L2GH*qdiv2                ! NN09 (34)
        Phi3 = Phi1 + e3c*L2GH*qdiv2                ! NN09 (35)
        Phi4 = Phi1 - e4c*L2GH*qdiv2                ! NN09 (36)
        Phi5 =        e5c*L2GM*qdiv2                ! NN09 (37)
        D_25 = max( 1.d-20, Phi2*Phi4 + Phi3*Phi5 ) ! NN09 (31)

        ! Compute stability functions
        if (q2 < q22) then
           Sm = qdiv*Sm2
           Sh = qdiv*Sh2
        else
           Sm = q2*A1*( Phi3 - 3.*C1*Phi4 )/D_25 ! NN09 (27)
           Sh = q2*A2*( Phi2 + 3.*C1*Phi5 )/D_25 ! NN09 (28)
        end if

        !
        itau(i,j,k) = q(i,j,k)/L(i,j,k)
!        Lq          = L(i,j,k)*q(i,j,k)
!        itau(i,j,k) = sqrt(2.*tke(i,j,k))/L(i,j,k)
        Lq          = L(i,j,k)*sqrt(2.*tke(i,j,k))

        ! Compute counter-gradient fluxes of conserved variables
        if ( MYNN_LEVEL == 2 ) then
           ! No counter-gradient flux
           whl_cg = 0.
           wqt_cg = 0.
           wb_cg  = 0.           
        else
           ! Compute thermodyanamic (co-)variances from level-2.5 closure
           hl2_25  = qdiv*B2*L2*Sh*dhldz**2.
           qt2_25  = qdiv*B2*L2*Sh*dqtdz**2.
           hlqt_25 = qdiv*B2*L2*Sh*dhldz*dqtdz
           
           ! Compute buoyancy (co-)variances
           hlthv_25 = Beta_hl(i,j,k)*hl2_25  + Beta_qt(i,j,k)*hlqt_25
           qtthv_25 = Beta_hl(i,j,k)*hlqt_25 + Beta_qt(i,j,k)*qt2_25
           thv2_25  = max(0.d0, Beta_hl(i,j,k)*hlthv_25 + Beta_qt(i,j,k)*qtthv_25)

           hlthv = Beta_hl(i,j,k)*hl2(i,j,k)  + Beta_qt(i,j,k)*hlqt(i,j,k)
           qtthv = Beta_hl(i,j,k)*hlqt(i,j,k) + Beta_qt(i,j,k)*qt2(i,j,k)
           thv2  = max(0.d0, Beta_hl(i,j,k)*hlthv + Beta_qt(i,j,k)*qtthv)

           ! Limit q2 so that L/q is less than 1/N for N2 > 0 (NN09 Section 2.7)
           if ( q2/L2 < -GH ) then
              q2 = -L2GH ! NN09 (56)
           end if

           !
           Phi2 = q2 - e2c*L2GH*qdiv2                                                  ! NN09 (34)
           D_p  = max( 1.d-20, Phi2*( Phi4 - Phi1 + q2 ) + Phi5*( Phi3 - Phi1 + q2 ) ) ! NN09 (32)
           EH   = qdiv*eHc*( Phi2 + Phi5 )/D_p                                         ! NN09 (48)

           ! Compute counter-gradient hl-flux
           hlthv_p = hlthv - hlthv_25
!!$           if (hlthv_25 >= 0.) then
!!$              hlthv_p = max( 0.d0, hlthv - hlthv_25 )
!!$           else
!!$              hlthv_p = min( 0.d0, hlthv - hlthv_25 )
!!$           end if
           whl_cg = Lq*EH*goth00*hlthv_p

           ! Compute counter-gradient qt-flux
           qtthv_p = qtthv - qtthv_25
!!$           if (qtthv_25 >= 0.) then
!!$              qtthv_p = max( 0.d0, qtthv - qtthv_25 )
!!$           else
!!$              qtthv_p = min( 0.d0, qtthv - qtthv_25 )
!!$           end if
           wqt_cg = Lq*EH*goth00*qtthv_p

           ! Restrict anisotropy by restricting buoyancy variance (NN09 Section 2.7)
           thv2_p = max( 0.d0, Beta_hl(i,j,k)*hlthv_p + Beta_qt(i,j,k)*qtthv_p )
           wden = ( 1. - C3 )*goth002*L2*qdiv2*( e4c*Phi2 - e3c*Phi5 )
           if (wden /= 0.) then
              Cw_25   = Phi1*( Phi2 + 3.*C1*Phi5 )/(3.*D_25) ! NN09 (59)
              Cw_low  = q2*( 0.12 - Cw_25 )*D_p/wden
              Cw_high = q2*( 0.76 - Cw_25 )*D_p/wden

              if (wden > 0.) then
                 thv2_p = min( max( thv2_p, Cw_low ), Cw_high )
              else
                 thv2_p = max( min( thv2_p, Cw_low ), Cw_high )
              end if
           end if

           ! Compute counter-gradient buoyancy flux
           wb_cg  = Lq*EH*goth002*thv2_p

           !
           EM = qdiv*eMc*( Phi3 - Phi4 )/(D_p*L2GH) ! NN09 (47)

           !
           Sm = Sm + EM*L2*goth002*thv2_p
        end if

        ! Compute turbulent diffusivities
        Km(i,j,k) = Lq*Sm
        Kh(i,j,k) = Lq*Sh

        !
        if (DOMF /= 0.) then
           whl = -Kh(i,j,k)*dhldz + whl_cg + whl_mf(i,j,k)
           wqt = -Kh(i,j,k)*dqtdz + wqt_cg + wqt_mf(i,j,k)
        else
           whl = -Kh(i,j,k)*dhldz + whl_cg
           wqt = -Kh(i,j,k)*dqtdz + wqt_cg
        end if
           
        ! Compute counter-gradient fluxes of GEOS variables
        ifac    = (zle(i,j,k) - zlo(i,j,k+1))*idzlo
        ac_half = ac(i,j,k+1) + ifac*(ac(i,j,k) - ac(i,j,k+1))
        wql     = ac_half*( A(i,j,k)*( -Kh(i,j,k)*dqtdz + wqt_cg ) - B(i,j,k)*( -Kh(i,j,k)*dhldz + whl_cg ) )

!        wql_cg(i,j,k) = 0.
        wql_cg(i,j,k) = wql + Kh(i,j,k)*dqldz
        ws_cg(i,j,k)  = MAPL_CP*whl_cg + MAPL_ALHL*wql_cg(i,j,k)
        wqv_cg(i,j,k) = wqt_cg - wql_cg(i,j,k)

        ! Compute budget terms for second-order moments
        tket_M(i,j,k) = Km(i,j,k)*S2(i,j,k)
        if (DOMF /= 0.) then
           tket_B(i,j,k) = -Kh(i,j,k)*N2(i,j,k) + wb_cg + goth00*wthv_mf(i,j,k)
        else
           tket_B(i,j,k) = -Kh(i,j,k)*N2(i,j,k) + wb_cg
        end if
        hl2t_M(i,j,k)  = -2.*whl*dhldz
        qt2t_M(i,j,k)  = -2.*wqt*dqtdz
        hlqtt_M(i,j,k) = -whl*dqtdz - wqt*dhldz

        ! Start test
        tau     = L(i,j,k)/q(i,j,k)

        w2      = q(i,j,k)**2./3. + 2.*A1*tau*( -Km(i,j,k)*S2(i,j,k) ) &
                                  + 4.*A1*( 1. - C2 )*tau*( -Kh(i,j,k)*N2(i,j,k) + wb_cg)
        whl_test = 3.*A2*tau*( -w2*dhldz     + ( 1. - C3 )*goth00*(  hlthv_25 + hlthv_p ) )
        wqt_test = 3.*A2*tau*( -w2*dqtdz     + ( 1. - C3 )*goth00*(  qtthv_25 + qtthv_p ) )
        wb_test  = 3.*A2*tau*( -w2*N2(i,j,k) + ( 1. - C3 )*goth002*( thv2_25  + thv2_p ) )
        Kh_test  = -( -Kh(i,j,k)*N2(i,j,k) + wb_cg )/N2(i,j,k)

        write(*,*) real(tke(i,j,k),4), real(w2/q2,4), real(qdiv,4), real(ac_half,4)

        ! End test
     end do
     end do
  end do

  ! Compute TKE diffusivity
  do k = 2,LM-1

     km1 = k - 1
     do j = 1,JM
     do i = 1,IM
        K_tke(i,j,k) = 1.5*( Km(i,j,k) + Km(i,j,km1) )
        K_tpe(i,j,k) = 0.5*( Km(i,j,k) + Km(i,j,km1) )
     end do
     end do
  end do
  K_tke(:,:,1)  = K_tke(:,:,2)
  K_tpe(:,:,1)  = K_tpe(:,:,2)
  K_tke(:,:,LM) = K_tke(:,:,LM-1)
  K_tpe(:,:,LM) = K_tpe(:,:,LM-1)

end subroutine run_mynn

!
! 
!
subroutine mynn_length(IM, JM, LM, wb_surf, zle, zlo, q, N2, LMO, L)

  use MAPL_ConstantsMod, only: MAPL_KARMAN

  double precision, parameter :: alpha1 = 0.23
  double precision, parameter :: alpha2 = 1.
  double precision, parameter :: alpha3 = 5.
  double precision, parameter :: alpha4 = 100.

  double precision, parameter :: qmin = 0.
  double precision, parameter :: zmax = 1.

  double precision, parameter :: cns = 2.7

  integer, intent(in)                                  :: IM, JM, LM
  double precision, intent(in)                         :: LMO
  double precision, dimension(IM,JM), intent(in)       :: wb_surf
  real, dimension(IM,JM,LM), intent(in)                :: zlo
  real, dimension(IM,JM,0:LM), intent(in)              :: zle
  double precision, dimension(IM,JM,0:LM), intent(in)  :: q, N2
  double precision, dimension(IM,JM,0:LM), intent(out) :: L

  integer                            :: i, j, k, kk, kkm1
  double precision                   :: qdz, LS, LB, LF, N, zeta
  double precision, dimension(IM,JM) :: work1, work2, LT, qc

  ! Compute integrals for turbulence length scale
  work1(:,:) = 1.E-5
  work2(:,:) = 1.E-5
  do k = 1,LM
     kk = LM - k + 1

     kkm1 = kk - 1
     do j = 1,JM
     do i = 1,IM
        qdz = max( 0., ( 0.5*( q(i,j,kk) + q(i,j,kkm1) ) - qmin )*( zle(i,j,kkm1) - zle(i,j,kk) ) )

        work1(i,j) = work1(i,j) + zlo(i,j,kk)*qdz 
        work2(i,j) = work2(i,j) + qdz
     end do
     end do
  end do

  ! Compute turbulence length scale and CBL vertical velocity scale
  do j = 1,JM
  do i = 1,IM
     LT(i,j) = alpha1*work1(i,j)/work2(i,j) ! NN09 (54)
     qc(i,j) = ( max( 0., wb_surf(i,j) )*LT(i,j) )**onethird
  end do
  end do

  ! Compute master length scale
  do k = 1,LM-1
     do j = 1,JM
     do i = 1,IM
        !
        zeta = zle(i,j,k)/LMO

        ! Compute surface length scale
!        ! NN09 (53)
!        if (zeta >= 1.) then
!           LS = MAPL_KARMAN*zle(i,j,k)/3.7
!        else if (zeta >= 0.) then
!           LS = MAPL_KARMAN*zle(i,j,k)/( 1. + cns*zeta )
!        else
!           LS = MAPL_KARMAN*zle(i,j,k)*( 1. - alpha4*zeta )**0.2
!        end if
        ! WRF version
        if (zeta > 0.) then
           LS = MAPL_KARMAN*zle(i,j,k)/( 1. + cns*min(zmax, zeta) )
        else
           LS = MAPL_KARMAN*zle(i,j,k)*( 1. - alpha4*zeta )**0.2
        end if
        
        ! Compute buoyancy length scale
        if (N2(i,j,k) > 0.) then
           N = sqrt(N2(i,j,k))
           
!           ! NN09 (55)
!           LF = alpha2*q(i,j,k)/N
!           if (zeta >= 0.) then
!              LB = LF
!           else
!              LB = ( alpha2*q(i,j,k) + alpha3*q(i,j,k)*sqrt(qc(i,j)/(LT(i,j)*N)) )/N
!           end if
           ! WRF version
           LB = alpha2*q(i,j,k)/N*( 1. + alpha3/alpha2*sqrt(qc(i,j)/(N*LT(i,j))) )
           LF = alpha2*q(i,j,k)/N
        else
           LB = 1.E+10
           LF = LB
        end if
        
        ! Harmonically average length scales
        L(i,j,k) = min( LF, LB/( LB/LS + LB/LT(i,j) + 1. ) ) ! NN09 (52)
     end do
     end do
  end do

end subroutine mynn_length

end module mynn
