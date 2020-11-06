       subroutine ocalbedo(rad,lam,aw,bw,wfac,sunz,ws,rod,ros)

!  Computes ocean surface albedo from solar zenith angle (sunz) 
!  and wind speed (ws, m/s).
!  Albedo is provided as direct (albd) and diffuse (albs).     

      integer, parameter  :: nlt=33
      integer lam(nlt)
      real aw(nlt),bw(nlt)
      real wfac(nlt)
      real rod(nlt),ros(nlt)

!  Derive surface reflectance as a function of sunz and ws
!  for OASIM Bands
      call sfcrfl(rad,lam,aw,bw,wfac,sunz,ws,rod,ros)

      return
      end
