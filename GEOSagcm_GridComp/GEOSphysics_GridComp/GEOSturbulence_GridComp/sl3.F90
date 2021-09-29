module sl3

use MAPL_ConstantsMod, only: mapl_grav, mapl_cp, mapl_alhl, mapl_rgas, mapl_rvap, mapl_p00, mapl_karman

implicit none

real, parameter :: th00 = 300.

real, parameter :: lvocp = mapl_alhl/mapl_cp
real, parameter :: kappa = mapl_rgas/MAPL_CP
real, parameter :: ep2   = mapl_rvap/mapl_rgas - 1.

real, parameter :: A1_sl3 = 0.92
real, parameter :: A2_sl3 = 0.74
real, parameter :: B1_sl3 = 16.6
real, parameter :: B2_sl3 = 10.1

! MYNN constants 
real, parameter ::             Pr     = 0.74
real, parameter ::             gamma1 = 0.235
real, parameter ::             B1     = 24.
real, parameter ::             B2     = 15.
real, parameter ::             C2     = 0.729
real, parameter ::             C4     = 0.
real, parameter ::             C5     = 0.2
double precision, parameter :: C3     = 0.34d0
!double precision, parameter :: C3     = 1.d0

double precision, parameter :: alpha1 = 0.23d0
double precision, parameter :: alpha2 = 1.d0
double precision, parameter :: alpha3 = 5.d0
double precision, parameter :: alpha4 = 100.d0

real, parameter :: A1 = B1*( 1. - 3.*gamma1 )/6.

double precision, parameter :: C1 = gamma1 - 1./(3.*A1*B1**(1./3.))

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

logical :: initialized = .false.

contains

!
! run_sl3
!
subroutine run_sl3(IM, JM, LM, &                                           
                   plo, z, zle, &                               
                   A_sl3, B_sl3, &                                         
                   u, v, s, qa, q, ql, qi, &                               
                   tkeshoc, qt2, hl2_sl3, hlqt_sl3, &                  
                   itau, km, kh, ws_explicit, wq_explicit, wql_explicit, & 
                   tket_M, tket_B, qt2t_M, hl2t_M, hlqtt_M)                

  integer, intent(in)                        :: IM, JM, LM
  real, dimension(IM,JM,LM), intent(in)      :: z, plo, u, v, qa, s, q, ql, qi, A_sl3, B_sl3
  real, dimension(IM,JM,0:LM), intent(in)    :: zle
  real, dimension(IM,JM,LM), intent(inout)   :: itau, tkeshoc, qt2, hl2_sl3, hlqt_sl3
  real, dimension(IM,JM,0:LM), intent(inout) :: km, kh, ws_explicit, wq_explicit, wql_explicit
  real, dimension(IM,JM,LM), intent(out)     :: tket_M, tket_B, qt2t_M, hl2t_M, hlqtt_M

  integer :: i, j, k, kp1, km1
  real :: idz_down, idz_up, dudz_down, dudz_up, dvdz_down, dvdz_up, &
          dhldz_down, dhldz_up, dqtdz_down, dqtdz_up, ifac, iexner, work, &
          goth00, whl_full, wqt_full, hlb, qtb, b2, dudz, dvdz, dhldz, dqtdz, &
          Beta_hl, Beta_qt, N2, S2 
  double precision :: LS, LT, L, sqrt_tke
  real, dimension(IM,JM,LM) :: hl, qt, km_full, kh_full, ws_full, wq_full, wql_full

  goth00 = mapl_grav/th00

  if ( .not. initialized ) then
     tkeshoc(:,:,:)  = 0.
     qt2(:,:,:)      = 0.
     hl2_sl3(:,:,:)  = 0.
     hlqt_sl3(:,:,:) = 0.

     tkeshoc(:,:,LM) = 0.1

     initialized = .true. 
  end if

  ! Initialize surface and TOA
  do j = 1,JM
  do i = 1,IM
     ! Surface 
     km(i,j,LM)           = 0.
     kh(i,j,LM)           = 0.
     ws_explicit(i,j,LM)  = 0.
     wq_explicit(i,j,LM)  = 0.
     wql_explicit(i,j,LM) = 0.

     ! TOA
     km(i,j,0)           = 0.
     kh(i,j,0)           = 0.
     ws_explicit(i,j,0)  = 0.
     wq_explicit(i,j,0)  = 0.
     wql_explicit(i,j,0) = 0.
  end do
  end do

  ! Compute conserved thermodynamic properties
  do k = 1,LM
     do j = 1,JM
     do i = 1,IM
        hl(i,j,k) = ( s(i,j,k) - mapl_alhl*ql(i,j,k) )/mapl_cp
        qt(i,j,k) = q(i,j,k) + ql(i,j,k) + qi(i,j,k)
     end do
     end do
  end do

  do k = 1,LM
     kp1 = k + 1
     km1 = k - 1
     do j = 1,JM
     do i = 1,IM
        ! Gradients at half levels
        if ( k >= 2 ) then
           idz_up = 1./( z(i,j,km1) - z(i,j,k) )

           dudz_up  = ( u(i,j,km1) - u(i,j,k) )*idz_up
           dvdz_up  = ( v(i,j,km1) - v(i,j,k) )*idz_up
           dhldz_up = ( hl(i,j,km1) - hl(i,j,k) )*idz_up
           dqtdz_up = ( qt(i,j,km1) - qt(i,j,k) )*idz_up 

           ifac = ( z(i,j,k) - zle(i,j,k) )/( zle(i,j,km1) - zle(i,j,k) )
        else
           dudz_up = 0.
           dvdz_up = 0.
           dhldz_up = 0.
           dqtdz_up = 0.

           ifac = 0.
        end if

        if ( k <= LM - 1 ) then
           idz_down = 1./( z(i,j,k) - z(i,j,kp1) )

           dudz_down  = ( u(i,j,k) - u(i,j,kp1) )*idz_down
           dvdz_down  = ( v(i,j,k) - v(i,j,kp1) )*idz_down
           dhldz_down = ( hl(i,j,k) - hl(i,j,kp1) )*idz_down
           dqtdz_down = ( qt(i,j,k) - qt(i,j,kp1) )*idz_down 
        else
           dudz_down  = dudz_up
           dvdz_down  = dvdz_up
           dhldz_down = dhldz_up
           dqtdz_down = dqtdz_up
        end if

        ! Gradients at full levels
        dudz  = dudz_down  + ifac*( dudz_up - dudz_down )        
        dvdz  = dvdz_down  + ifac*( dvdz_up - dvdz_down )        
        dhldz = dhldz_down + ifac*( dhldz_up - dhldz_down )
        dqtdz = dqtdz_down + ifac*( dqtdz_up - dqtdz_down )
     
        ! Moist-turbulence coefficients
        iexner = ( MAPL_P00/plo(i,j,k) )**kappa
        work   = lvocp*iexner - ( 1. + ep2 )*th00

        Beta_hl = iexner*( 1. - work*qa(i,j,k)*B_sl3(i,j,k) )
        Beta_qt = ep2*th00 + work*qa(i,j,k)*A_sl3(i,j,k)

        ! Turbulence parameters
        S2 = dudz**2. + dvdz**2.

        N2  = goth00*( Beta_hl*dhldz + Beta_qt*dqtdz )
        hlb = goth00*( Beta_hl*hl2_sl3(i,j,k)  + Beta_qt*hlqt_sl3(i,j,k) )
        qtb = goth00*( Beta_hl*hlqt_sl3(i,j,k) + Beta_qt*qt2(i,j,k) )

        b2 = goth00*( Beta_hl*hlb + Beta_qt*qtb )

        ! Length scale
        LT = 400.
        LS = mapl_karman*z(i,j,k)
        L  = LT/( 1. + LT/LS )

        sqrt_tke = sqrt(tkeshoc(i,j,k))

        ! Diffusivities
!        km_full(i,j,k) = A1_sl3*l*sqrt_tke
!        kh_full(i,j,k) = A2_sl3*l*sqrt_tke
        call mynn_l25(tkeshoc(i,j,k), L, S2, N2, &
                      km_full(i,j,k), kh_full(i,j,k))

        ! Explicit fluxes
        whl_full = 0.
        wqt_full = 0.
        
        wql_full(i,j,k) = 0. 
        ws_full(i,j,k)  = mapl_cp*whl_full + mapl_alhl*wql_full(i,j,k)
        wq_full(i,j,k)  = wqt_full - wql_full(i,j,k)

        ! Outputs variables
        itau(i,j,k)    = sqrt_tke/L
        tket_M(i,j,k)  = km_full(i,j,k)*S2
        tket_B(i,j,k)  = -kh_full(i,j,k)*N2
        qt2t_M(i,j,k)  = kh_full(i,j,k)*dqtdz**2.
        hl2t_M(i,j,k)  = kh_full(i,j,k)*dhldz**2.
        hlqtt_M(i,j,k) = kh_full(i,j,k)*dhldz*dqtdz        

        ! Debug output
        write(*,*) k, tkeshoc(i,j,k), sqrt(1.E+6*qt2(i,j,k)), qa(i,j,k)
     end do
     end do
  end do

  do k = 2,LM-1

     km1 = k - 1
     kp1 = k + 1
     do j = 1,JM
     do i = 1,IM
        ! Output variables
        km(i,j,k)           = 0.5*( km_full(i,j,k)  + km_full(i,j,kp1) )
        kh(i,j,k)           = 0.5*( kh_full(i,j,k)  + kh_full(i,j,kp1) )
        ws_explicit(i,j,k)  = 0.5*( ws_full(i,j,k)  + ws_full(i,j,kp1) )
        wq_explicit(i,j,k)  = 0.5*( wq_full(i,j,k)  + wq_full(i,j,kp1) )
        wql_explicit(i,j,k) = 0.5*( wql_full(i,j,k) + wql_full(i,j,kp1) )
     end do
     end do
  end do

end subroutine run_sl3

!
! mynn_l25
!
subroutine mynn_l25(tke, L, S2, N2, &
                    km, kh)

  implicit none

  ! Inputs/outputs
  real, intent(in)             :: tke, S2, N2 
  double precision, intent(in) :: L
  real, intent(out)            :: km, kh

  ! Local variables
  double precision :: Phi1, Phi2, Phi3, Phi4, Phi5, D_25, qdiv, qdiv2, &
                      q22, L2, q2, GH, GM, L2GH, L2GM, SM2, SH2, q, Lq, &
                      SM25, SH25

  !
  q = sqrt( max( 1.d-10, real(2.*tke,8) ) )

  ! Compute some intermediate quantities
  Lq   = L*q
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

  ! Set diffusivities
  km = Lq*SM25
  kh = Lq*SH25

end subroutine mynn_l25

!
! mynn_l2
!
subroutine mynn_l2(GM, GH, SM2, SH2)

  double precision, intent(in)  :: GM, GH
  double precision, intent(out) :: SM2, SH2

  double precision :: Ri, Rf

  Ri  = -GH/max( GM, 1.d-10 )
  Rf  = min( Rfc, Ri1*( Ri + Ri2 - sqrt( Ri**2.d0 - Ri3*Ri + Ri4 ) ) ) ! NN09 (A11)

  SH2 = SHc*( ( Rfc - Rf )/( 1.d0 - Rf ) )     ! NN09 (A4)
  SM2 = SMc*( ( Rf1 - Rf )/( Rf2 - Rf )  )*SH2 ! NN09 (A3)

end subroutine mynn_l2

end module sl3
