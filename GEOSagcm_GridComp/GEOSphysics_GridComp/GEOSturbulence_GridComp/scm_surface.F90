module scm_surface

use MAPL_ConstantsMod, only: mapl_grav, mapl_cp, mapl_karman, mapl_rdry, mapl_pi, mapl_p00
use MAPL_SatVaporMod, only: MAPL_EQsat

implicit none

real, parameter :: a = 15.
real, parameter :: b = 4.7

logical :: initialized = .false.

contains

subroutine surface(IM, JM, LM, &                            ! in
                   surf_type, T_surf, rh_surf, dTdt_surf, & ! in
                   dt, ple, &                               ! in 
                   s_surf, &                                ! inout
                   q_surf)                                  ! out

  integer, intent(in)                     :: IM, JM, LM, surf_type
  real, intent(in)                        :: dt, T_surf, rh_surf, dTdt_surf
  real, dimension(IM,JM,0:LM), intent(in) :: ple
  real, dimension(IM,JM), intent(inout)   :: s_surf
  real, dimension(IM,JM), intent(out)     :: q_surf

  integer :: i, j, k

  if ( .not. initialized ) then
     s_surf(:,:) = mapl_cp*T_surf
     
     initialized = .true.
  end if
  
  do j = 1,JM
  do i = 1,IM
     if ( surf_type == 1 ) then
        s_surf(i,j) = s_surf(i,j) + dt*mapl_cp*dTdt_surf
     end if

     q_surf(i,j) = rh_surf*MAPL_EQsat(s_surf(i,j)/mapl_cp, ple(i,j,LM)) 
  end do
  end do

end subroutine surface
                   

!
! surface_layer
!
subroutine surface_layer(IM, JM, LM, &
                         flux_type, z0, &
                         zi, s_surf, q_surf, &
                         zl, zle, ple, rhoe, u, v, T, qv, thv, &
                         sh, evap, zeta, & 
                         u_star, cm, ct)

  integer, intent(in)                     :: IM, JM, LM, flux_type
  real, intent(in)                        :: z0
  real, dimension(IM,JM), intent(in)      :: zi, s_surf, q_surf
  real, dimension(IM,JM,0:LM), intent(in) :: zle, ple, rhoe
  real, dimension(IM,JM,LM), intent(in)   :: zl, u, v, T, qv, thv
  real, dimension(IM,JM), intent(inout)   :: sh, evap, zeta
  real, dimension(IM,JM), intent(out)     :: u_star, cm, ct

  integer :: i, j, k, iter, niter
  real :: work1, work2, Ri_bulk, foo, bar, Dfoo, Dbar, f, Df
  
  real, dimension(IM,JM) :: w_star, bw_surf, wspd

  if ( flux_type == 2 ) then
     niter = 1
  else
     niter = 20

     zeta(:,:) = -0.01
     
     if ( flux_type == 1 ) then
        bw_surf(:,:) = (mapl_grav/thv(:,:,LM)) * ( sh(:,:)/mapl_cp + evap(:,:) )/rhoe(:,:,LM)
     else
        bw_surf(:,:) = 0.01
     end if
  end if

  !
  do iter = 1,niter
     do j = 1,JM
     do i = 1,IM
        ! Compute w_star
        if ( bw_surf(i,j) > 0. ) then
           w_star(i,j) = ( bw_surf(i,j)*max( zi(i,j), zle(i,j,LM-1) ) )**(1./3.)
        else
           w_star(i,j) = 0.
        end if

        ! Compute windspeed
        wspd(i,j) = sqrt( u(i,j,LM)**2. + v(i,j,LM)**2. + w_star(i,j)**2. )

        ! Approximate zeta using Newton's method
        work1 = ( zl(i,j,LM) + z0 )/z0
        work2 = zl(i,j,LM) + z0

        foo = log(work1) - Psi_H(work2/zeta(i,j)) + Psi_H(z0/zeta(i,j))
        bar = log(work1) - Psi_M(work2/zeta(i,j)) + Psi_M(z0/zeta(i,j))

        ! Approximate zeta using Newton's method
        if ( flux_type /= 2 ) then
           Dfoo = ( work2*DPsi_H(work2/zeta(i,j)) - z0*DPsi_H(z0/zeta(i,j)) )/zeta(i,j)**2.
           Dbar = ( work2*DPsi_M(work2/zeta(i,j)) - z0*DPsi_M(z0/zeta(i,j)) )/zeta(i,j)**2.

           Ri_bulk = zl(i,j,LM)*bw_surf(i,j)/wspd(i,j)**3.

           f  = mapl_karman**2.*zeta(i,j)*foo/bar**2. - Ri_bulk
           Df = mapl_karman**2.*( ( foo + zeta(i,j)*Dfoo )/bar**2. - 2.*zeta(i,j)*foo*Dbar/bar**3. )

           zeta(i,j) = zeta(i,j) - f/Df
        end if

        ! Compute surface thermodynamic fluxes
        if ( flux_type /= 1 ) then 
           ct(i,j) = mapl_karman**2./( foo*bar )

           sh(i,j)   = -ct(i,j)*rhoe(i,j,LM)*( mapl_cp*T(i,j,LM) + mapl_grav*zl(i,j,LM) - s_surf(i,j) )
           evap(i,j) = -ct(i,j)*rhoe(i,j,LM)*( qv(i,j,LM) - q_surf(i,j) )

           bw_surf(i,j) = (mapl_grav/thv(i,j,LM)) * ( sh(i,j)/mapl_cp + thv(i,j,LM)*evap(i,j) )/rhoe(i,j,LM) 
        else
           ct(i,j) = 0.
        end if

        ! Compute momentum exchange coefficient
        cm(i,j) = mapl_karman**2./bar**2.

        ! Compute u_star
        u_star(i,j) = sqrt( cm(i,j) )*wspd(i,j)
     end do
     end do
  end do
 
end subroutine surface_layer

!
! Phi_M
!
real function Phi_M(zeta)

  real, intent(in) :: zeta

  if ( zeta < 0. ) then
     Phi_M = ( 1. - a*zeta )**(-0.25)
  else
     Phi_M = 1. + b*zeta
  end if

end function Phi_M

!
! Phi_H
!
real function Phi_H(zeta)

  real, intent(in) :: zeta

  if ( zeta < 0. ) then
     Phi_H = ( 1. - a*zeta )**(-0.5)
  else
     Phi_H = 1. + b*zeta
  end if

end function Phi_H

!
! Psi_M
!
real function Psi_M(zeta)

  real, intent(in) :: zeta

  real :: x

  if ( zeta < 0. ) then
     x     = ( 1. - a*zeta )**0.25
     Psi_M = 0.5*mapl_pi - 2.*atan(x) + log( 0.125*( 1. + x )**2.*( 1. + x**2. ) )
  else
     Psi_M = -b*zeta
  end if

end function Psi_M

!
! Psi_H
!
real function Psi_H(zeta)

  real, intent(in) :: zeta

  real :: x

  if ( zeta < 0. ) then
     x     = ( 1. - a*zeta )**0.25
     Psi_H = 2.*log( 0.5*( 1. + x**2. ) )
  else
     Psi_H = -b*zeta
  end if

end function Psi_H

!
! DPsi_M
!
real function DPsi_M(zeta)

  real, intent(in) :: zeta

  DPsi_M = ( 1. - Phi_M(zeta) )/zeta

end function DPsi_M

!
! DPsi_H
!
real function DPsi_H(zeta)

  real, intent(in) :: zeta

  DPsi_H = ( 1. - Phi_H(zeta) )/zeta

end function DPsi_H

end module scm_surface
