module sl3

implicit none

contains

subroutine run_sl3(IM, JM, LM, &
                   zle, zlo, ple, plo, &
                   u, v, T, qv, ql, tke, hl2, qt2, hlqt, & ! State variables
                   ustar, ac, whl_mf, wqt_mf, wthv_mf, &
                   Kh, Km, Ke, itau, &
                   ws_cg, wqv_cg, wql_cg, &
                   tket_M, tket_B, hl2t_M, qt2t_M, hlqtt_M)

  use MAPL_SatVaporMod,  only: MAPL_EQsat
  use MAPL_ConstantsMod, only: MAPL_GRAV, MAPL_CP, MAPL_ALHL, MAPL_RGAS, MAPL_RVAP

  implicit none

  ! Constant parameters
  real, parameter :: vonk   = 0.4
  real, parameter :: tau_bl = 400.
  real, parameter :: zi     = 1000.
  real, parameter :: gocp   = MAPL_GRAV/MAPL_CP
  real, parameter :: lvocp  = MAPL_ALHL/MAPL_CP
  real, parameter :: kappa  = MAPL_RGAS/MAPL_CP
  real, parameter :: Pr     = 0.74
  real, parameter :: th00   = 300.
  real, parameter :: Ce     = 0.16
  real, parameter :: beta   = MAPL_GRAV/th00
  real, parameter :: ep2    = MAPL_RVAP/MAPL_RGAS - 1.

  real, parameter :: A1 = 1.18 
  real, parameter :: A2 = 0.665
  real, parameter :: B1 = 24.
  real, parameter :: B2 = 15.
  real, parameter :: C3 = 0.294

  ! Inputs
  integer, intent(in) :: IM, JM, LM

  real, dimension(IM,JM,0:LM), intent(in) :: zle  ! Height at grid cell interfaces [m]
  real, dimension(IM,JM,LM), intent(in)   :: zlo  ! Height at grid cell centers [m]
  real, dimension(IM,JM,0:LM), intent(in) :: ple  ! Pressure of grid cell interfaces [Pa]
  real, dimension(IM,JM,LM), intent(in)   :: plo  ! Pressure at grid cell centers [Pa]
  real, dimension(IM,JM,LM), intent(in)   :: u    ! Mean east-west wind speed [ms-1]
  real, dimension(IM,JM,LM), intent(in)   :: v    ! Mean north-south wind speed [ms-1]
  real, dimension(IM,JM,LM), intent(in)   :: T    ! Mean temperature [K]
  real, dimension(IM,JM,LM), intent(in)   :: qv   ! Mean water vapor specific humidity [kgkg-1] 
  real, dimension(IM,JM,LM), intent(in)   :: ql   ! Mean liquid water specific humidity [kgkg-1]
  real, dimension(IM,JM,0:LM), intent(inout) :: tke  ! Turbulent kinetic energy [m+2s-2]
  real, dimension(IM,JM,0:LM), intent(in) :: hl2  ! Variance of liquid water static energy [K+2]
  real, dimension(IM,JM,0:LM), intent(in) :: qt2  ! Variance of total water specific humidity [kg+2kg-2]
  real, dimension(IM,JM,0:LM), intent(in) :: hlqt ! Covariance of hl and qt [Kkgkg-1]
  real, dimension(IM,JM), intent(in)      :: ustar
  real, dimension(IM,JM,LM), intent(in)   :: ac   ! Subgrid cloud fraction [ ]
  real, dimension(IM,JM,0:LM), intent(in) :: whl_mf   
  real, dimension(IM,JM,0:LM), intent(in) :: wqt_mf  
  real, dimension(IM,JM,0:LM), intent(in) :: wthv_mf  

  ! Outputs
  real, dimension(IM,JM,0:LM), intent(out) :: Kh      ! Eddy conductivity [m+2s-1]
  real, dimension(IM,JM,0:LM), intent(out) :: Km      ! Eddy diffusivity [m+2s-1]
  real, dimension(IM,JM,LM), intent(out)   :: Ke      ! Diffusivity coefficient for TKE [m+2s-1]
  real, dimension(IM,JM,0:LM), intent(out) :: itau    ! Inverse of TKE dissipation timescale [s-1]
  real, dimension(IM,JM,0:LM), intent(out) :: ws_cg
  real, dimension(IM,JM,0:LM), intent(out) :: wqv_cg
  real, dimension(IM,JM,0:LM), intent(out) :: wql_cg
  real, dimension(IM,JM,0:LM), intent(out) :: tket_M  ! TKE mean-gradient production rate [m+2s-3] 
  real, dimension(IM,JM,0:LM), intent(out) :: tket_B  ! TKE buoyancy production rate [m+2s-3]
  real, dimension(IM,JM,0:LM), intent(out) :: hl2t_M  ! hl2 mean-gradient production rate [K+2s-1]
  real, dimension(IM,JM,0:LM), intent(out) :: qt2t_M  ! qt2 mean-gradient prodcution rate [kg+2kg-2s-1]
  real, dimension(IM,JM,0:LM), intent(out) :: hlqtt_M ! hlqt mean-gradient production rate [Kkgkg-1]

  ! Local arrays
  real, dimension(IM,JM,LM) :: hl, qt, thl

  ! Local scalars
  integer :: i, j, k, kp1, km1
  real    :: l, l1, l2, l3, l23, dzlo, dhldz, dqtdz, dqldz, ac_half, A, B, T_half, ql_half, &
             whl, wqt, ifac, N2, S2, exner, q, iexner, Tl, s, qs, dqs, whl_cg, wqt_cg, wql, &
             hlthv, qtthv, thv2, wb_cg, w2, tau
  real    :: wrk1, wrk2, wrk3

  ! Test
  tke(:,:,LM-1) = max( 1.E-3, tke(:,:,LM-1) )

  ! Compute conserved thermodyanamic properties
  do k = 1,LM
     do j = 1,JM
     do i = 1,IM
        hl(i,j,k) = T(i,j,k) + gocp*zlo(i,j,k) - lvocp*ql(i,j,k)
        qt(i,j,k) = qv(i,j,k) + ql(i,j,k)

        iexner     = (1E+5/plo(i,j,k))**kappa
        thl(i,j,k) = iexner*( T(i,j,k) - lvocp*ql(i,j,k) )
     end do
     end do
  end do
  
  !
  do k = 1,LM-1

     kp1 = k + 1
     do j = 1,JM
     do i = 1,IM
        ! Compute gradients
        dzlo   = zlo(i,j,k) - zlo(i,j,kp1)
        dhldz  = (hl(i,j,k) - hl(i,j,kp1))/dzlo
        dqtdz  = (qt(i,j,k) - qt(i,j,kp1))/dzlo
        dqldz  = (ql(i,j,k) - ql(i,j,kp1))/dzlo

        ! Compute useful quantities
        exner    = (ple(i,j,k)*1.E-5)**kappa
        ifac     = (zle(i,j,k) - zlo(i,j,kp1))/dzlo
        iexner   = 1./exner

        ! Interpolation stuff to half levels
        ac_half = ac(i,j,kp1) + ifac*(ac(i,j,k) - ac(i,j,kp1))
        T_half  = T(i,j,kp1)  + ifac*(T(i,j,k)  - T(i,j,kp1))
        ql_half = ql(i,j,kp1) + ifac*(ql(i,j,k) - ql(i,j,kp1))

        !
        Tl = T_half - lvocp*ql_half
        qs = MAPL_EQsat(Tl, ple(i,j,k), dqs)
        A  = 1./( 1. + lvocp*dqs )
        B  = A*exner*dqs

        wrk1 = lvocp*iexner - (1. + ep2)*th00
        wrk2 = 1. - ac_half*B*wrk1 
        wrk3 = ep2*th00 + ac_half*A*wrk1 

        ! Compute buoyancy frequency
        N2 = beta*( wrk2*dhldz*iexner + wrk3*dqtdz )

        ! Compute shear
        S2 = ((u(i,j,k) - u(i,j,kp1))/dzlo)**2. + ((v(i,j,k) - v(i,j,kp1))/dzlo)**2.

        ! Compute turbulent length and time scales
        q  = sqrt(2.*tke(i,j,k))
        l1 = vonk*zle(i,j,k)
        l2 = 400.
        l23 = 1./( 1./l1 + 1./l2  )
        if (q > 0.) then
           if (N2 > 0.) then
              l3 = q/sqrt(N2)
              l  = 1./( 1./l23 + 1./l3 )
           else
              l = l23 
           end if
           tau         = l/q
           itau(i,j,k) = q/l
        else
           if (N2 > 0.) then
              tau         = 1./sqrt(N2)
              itau(i,j,k) = sqrt(N2)
           else
              tau         = 0.
              itau(i,j,k) = 0.
           end if
        end if

        ! Temporary simplification (isotropic approximation)
        w2 = (2./3.)*tke(i,j,k)

        ! Compute diffusivities
        Kh(i,j,k) = 3.*A2*tau*w2
        Km(i,j,k) = 3.*A1*tau*w2

        ! Compute countergradient fluxes of conserved variables
        hlthv = wrk2*hl2(i,j,k)*iexner  + wrk3*hlqt(i,j,k)
        qtthv = wrk2*hlqt(i,j,k)*iexner + wrk3*qt2(i,j,k)
        thv2  = wrk2*hlthv*iexner       + wrk3*qtthv

!        whl_cg = 3.*A2*(1. - C3)*tau*beta*hlthv
!        wqt_cg = 3.*A2*(1. - C3)*tau*beta*qtthv
!        wb_cg  = 3.*A2*(1. - C3)*tau*beta**2.*thv2
        whl_cg = 0.
        wqt_cg = 0.
        wb_cg  = 0.

        ! Compute liquid water flux
        wql = ac_half*( A*( -Kh(i,j,k)*dqtdz + wqt_cg ) - B*( -Kh(i,j,k)*dhldz + whl_cg  )*iexner )

        ! Compute countergradient fluxes of prognostic thermodynamic state variables
!        wql_cg(i,j,k) = wql + Kh(i,j,k)*dqldz
!        ws_cg(i,j,k)  = mapl_cp*whl_cg + mapl_alhl*wql_cg(i,j,k)
!        wqv_cg(i,j,k) = wqt_cg - wql_cg(i,j,k)
        wql_cg(i,j,k) = 0.
        ws_cg(i,j,k)  = 0.
        wqv_cg(i,j,k) = 0.

        ! Compute fluxes
        whl = -Kh(i,j,k)*dhldz + whl_cg + whl_mf(i,j,k)
        wqt = -Kh(i,j,k)*dqtdz + wqt_cg + wqt_mf(i,j,k)

        ! Compute budget terms for second-order moments
        tket_M(i,j,k)  = Km(i,j,k)*S2
        tket_B(i,j,k)  = -Kh(i,j,k)*N2 + wb_cg + beta*wthv_mf(i,j,k)
        hl2t_M(i,j,k)  = -2.*whl*dhldz
        qt2t_M(i,j,k)  = -2.*wqt*dqtdz
        hlqtt_M(i,j,k) = -whl*dqtdz - wqt*dhldz

        write(*,*) tke(i,j,k), -Kh(i,j,k)*N2, beta*wthv_mf(i,j,k)
     end do
     end do
  end do

  ! 
  do k = 2,LM-1

     km1 = k - 1
     do j = 1,JM
     do i = 1,IM
        Ke(i,j,k) = 0.5*( Km(i,j,k) + Km(i,j,km1) ) 
     end do
     end do
  end do
  Ke(:,:,1)  = Ke(:,:,2)
  Ke(:,:,LM) = Ke(:,:,LM-1)

end subroutine run_sl3

end module sl3
