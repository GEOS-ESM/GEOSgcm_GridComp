module ncar_ocean_fluxes_mod

implicit none
private

public :: ncar_ocean_fluxes

real, parameter :: GRAV   = 9.80
real, parameter :: VONKARM = 0.40

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following Large and Yeager (used in NCAR models)           !
! Coded by Mike Winton (Michael.Winton@noaa.gov) in 2004
!
! A bug was found by Laurent Brodeau (brodeau@gmail.com) in 2007.
! Stephen.Griffies@noaa.gov updated the code with the bug fix. 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
subroutine ncar_ocean_fluxes (u_del, t, ts, q, qs, z, avail, &
                              cd, ch, ce, ustar, bstar       )
real   , intent(in)   , dimension(:) :: u_del, t, ts, q, qs, z
logical, intent(in)   , dimension(:) :: avail
real   , intent(inout), dimension(:) :: cd, ch, ce, ustar, bstar

  real :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real :: cd_rt                                ! full drag coefficients @ z
  real :: zeta, x2, x, psi_m, psi_h            ! stability parameters
  real :: u, u10, tv, tstar, qstar, z0, xx, stab
  integer, parameter :: n_itts = 2
  integer               i, j


  do i=1,size(u_del)
    if (avail(i)) then
      tv = t(i)*(1+0.608*q(i));
      u = max(u_del(i), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
      u10 = u;                                                ! first guess 10m wind
    
      cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
      cd_n10_rt = sqrt(cd_n10);
      ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
      stab = 0.5 + sign(0.5,t(i)-ts(i))
      ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c
  
      cd(i) = cd_n10;                                         ! first guess for exchange coeff's at z
      ch(i) = ch_n10;
      ce(i) = ce_n10;
      do j=1,n_itts                                           ! Monin-Obukhov iteration
        cd_rt = sqrt(cd(i));
        ustar(i) = cd_rt*u;                                   ! L-Y eqn. 7a
        tstar    = (ch(i)/cd_rt)*(t(i)-ts(i));                ! L-Y eqn. 7b
        qstar    = (ce(i)/cd_rt)*(q(i)-qs(i));                ! L-Y eqn. 7c
        bstar(i) = grav*(tstar/tv+qstar/(q(i)+1/0.608));
        zeta     = vonkarm*bstar(i)*z(i)/(ustar(i)*ustar(i)); ! L-Y eqn. 8a
        zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
        x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
        x2 = max(x2, 1.0);                                    ! undocumented NCAR
        x = sqrt(x2);
    
        if (zeta > 0) then
          psi_m = -5*zeta;                                    ! L-Y eqn. 8c
          psi_h = -5*zeta;                                    ! L-Y eqn. 8c
        else
          psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
          psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
        end if
    
        u10 = u/(1+cd_n10_rt*(log(z(i)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9
        cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
        cd_n10_rt = sqrt(cd_n10);
        ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
        stab = 0.5 + sign(0.5,zeta)
        ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
        z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic
    
        xx = (log(z(i)/10)-psi_m)/vonkarm;
        cd(i) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
        xx = (log(z(i)/10)-psi_h)/vonkarm;
!!$        ch(i) = ch_n10/(1+ch_n10*xx/cd_n10_rt)**2;                !       b (bug)
!!$        ce(i) = ce_n10/(1+ce_n10*xx/cd_n10_rt)**2;                !       c (bug)
        ch(i) = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd(i)/cd_n10) ! 10b (corrected code aug2007)
        ce(i) = ce_n10/(1+ce_n10*xx/cd_n10_rt)*sqrt(cd(i)/cd_n10) ! 10c (corrected code aug2007)
      end do
    end if
  end do

end subroutine ncar_ocean_fluxes



end module ncar_ocean_fluxes_mod
