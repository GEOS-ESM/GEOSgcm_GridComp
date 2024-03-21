      subroutine sfcrfl(rad,lam,aw,bw,wfac,theta,wspd,rod,ros)

!  Computes surface reflectance for direct (rod) and diffuse (ros)
!  components separately, as a function of theta, wind speed or
!  stress.
!  Includes spectral dependence of foam reflectance derived from Frouin
!  et al., 1996 (JGR)

      integer, parameter :: nlt=33
      save ifst,rn,roair
      integer lam(nlt)
      real aw(nlt),bw(nlt)
      real rod(nlt),ros(nlt)
      real wfac(nlt)
      data a0,a1,a2,a3 /0.9976,0.2194,5.554E-2,6.7E-3/
      data b0,b1,b2,b3 /5.026,-0.01138,9.552E-6,-2.698E-9/
      data ifst /0/

      if (ifst .eq. 0)then
       rn = 1.341    !index of refraction of pure seawater
       roair = 1.2E3     !density of air g/m3
!   This section computed in setlte.F90
!       do nl = 1,nlt
!        rlam = float(lam(nl))
!        if (lam(nl) .lt. 900)then
!         t = exp(-(aw(nl)+0.5*bw(nl)))
!         tlog = alog(1.0E-36+t)
!         fac = a0 + a1*tlog + a2*tlog*tlog + a3*tlog*tlog*tlog
!         wfac(nl) = min(fac,1.0)
!         wfac(nl) = max(fac,0.0)
!        else
!         fac = b0 + b1*rlam + b2*rlam*rlam + b3*rlam*rlam*rlam
!         wfac(nl) = max(fac,0.0)
!        endif
!       enddo
       ifst = 1
      endif

!  Foam and diffuse reflectance
      if (wspd .gt. 4.0)then
       if (wspd .le. 7.0)then
        cn = 6.2E-4 + 1.56E-3/wspd
        rof = roair*cn*2.2E-5*wspd*wspd - 4.0E-4
       else
        cn = 0.49E-3 + 0.065E-3*wspd
        rof = (roair*cn*4.5E-5 - 4.0E-5)*wspd*wspd
       endif
       rosps = 0.057
      else
       rof = 0.0
       rosps = 0.066
      endif

!  Direct
!   Fresnel reflectance for theta < 40, wspd < 2 m/s
      if (theta .lt. 40.0 .or. wspd .lt. 2.0)then
       if (theta .eq. 0.0)then
        rospd = 0.0211
       else
        rtheta = theta/rad
        sintr = sin(rtheta)/rn
        rthetar = asin(sintr)
        rmin = rtheta - rthetar
        rpls = rtheta + rthetar
        sinrmin = sin(rmin)
        sinrpls = sin(rpls)
        tanrmin = tan(rmin)
        tanrpls = tan(rpls)
        sinp = (sinrmin*sinrmin)/(sinrpls*sinrpls)
        tanp = (tanrmin*tanrmin)/(tanrpls*tanrpls)
        rospd = 0.5*(sinp + tanp)
       endif
      else
!   Empirical fit otherwise
       a = 0.0253
       b = -7.14E-4*wspd + 0.0618
       rospd = a*exp(b*(theta-40.0))
      endif

!  Reflectance totals
      do nl = 1,nlt
       rod(nl) = rospd + rof*wfac(nl)
       ros(nl) = rosps + rof*wfac(nl)
      enddo

      return
      end

