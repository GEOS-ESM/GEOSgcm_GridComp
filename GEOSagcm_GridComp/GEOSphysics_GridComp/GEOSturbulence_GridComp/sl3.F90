module sl3

use MAPL_ConstantsMod, only: mapl_grav, mapl_cp, mapl_alhl, mapl_rgas, mapl_rvap, mapl_p00, mapl_karman

implicit none

real, parameter :: th00 = 300.

real, parameter :: lvocp = mapl_alhl/mapl_cp
real, parameter :: kappa = mapl_rgas/MAPL_CP
real, parameter :: ep2   = mapl_rvap/mapl_rgas - 1.

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
                   tkeshear, tkebuoy, qt2t_M, hl2t_M, hlqtt_M)                

  integer, intent(in)                        :: IM, JM, LM
  real, dimension(IM,JM,LM), intent(in)      :: z, plo, u, v, qa, s, q, ql, qi, A_sl3, B_sl3
  real, dimension(IM,JM,0:LM), intent(in)    :: zle
  real, dimension(IM,JM,LM), intent(inout)   :: itau, tkeshoc, qt2, hl2_sl3, hlqt_sl3
  real, dimension(IM,JM,0:LM), intent(inout) :: km, kh, ws_explicit, wq_explicit, wql_explicit
  real, dimension(IM,JM,LM), intent(out)     :: tkeshear, tkebuoy, qt2t_M, hl2t_M, hlqtt_M

  integer :: i, j, k, kp1, km1
  real :: idz_down, idz_up, dudz_down, dudz_up, dvdz_down, dvdz_up, &
          dhldz_down, dhldz_up, dqtdz_down, dqtdz_up, ifac, iexner, work, &
          goth00, whl_full, wqt_full, hlb, qtb, b2, dudz, dvdz, dhldz, dqtdz, &
          Beta_hl, Beta_qt, N2, S2 
  double precision :: ls, lt, l, sqrt_tke
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
           idz_down = 1./( z(i,j,km1) - z(i,j,k) )

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
        lt = 400.
        ls = mapl_karman*z(i,j,k)
        l  = lt/( 1. + lt/ls )

        sqrt_tke = sqrt(tkeshoc(i,j,k))

        ! Diffusivities
        km_full(i,j,k) = l*sqrt_tke
        kh_full(i,j,k) = km(i,j,k)

        ! Explicit fluxes
        whl_full = 0.
        wqt_full = 0.
        
        wql_full(i,j,k) = 0. 
        ws_full(i,j,k)  = mapl_cp*whl_full + mapl_alhl*wql_full(i,j,k)
        wq_full(i,j,k) = wqt_full - wql_full(i,j,k)

        ! Outputs variables
        itau(i,j,k)     = sqrt_tke/l
        tkeshear(i,j,k) = km_full(i,j,k)*S2
        tkebuoy(i,j,k)  = -kh_full(i,j,k)*N2
        qt2t_M(i,j,k)   = kh_full(i,j,k)*dqtdz**2.
        hl2t_M(i,j,k)   = kh_full(i,j,k)*dhldz**2.
        hlqtt_M(i,j,k)  = kh_full(i,j,k)*dhldz*dqtdz        

        ! Debug output
        write(*,*) k, km_full(i,j,k), qa(i,j,k)
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

end module sl3
