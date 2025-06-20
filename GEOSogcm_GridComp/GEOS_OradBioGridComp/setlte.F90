      subroutine setlte(rad,pi2,lam,aw,bw,ac,bc,bpic,excdom,exdet,    &
                        WtoQ,wfac)

!  Sets parameters for ocean irradiance.

#include "definebio.h"
#include "comlte.h"
      integer lam(nlt)
      real aw(nlt),bw(nlt)
      real ac(nchl,nlt),bc(nchl,nlt)
      real bpic(nlt)
      real excdom(nlt),exdet(nlt)
      real WtoQ(nlt)
      real wfac(nlt)
      data a0,a1,a2,a3 /0.9976,0.2194,5.554E-2,6.7E-3/
      data b0,b1,b2,b3 /5.026,-0.01138,9.552E-6,-2.698E-9/

!  Degrees to radians conversion
      pi = dacos(-1.0D0)
      pi2 = pi*2.0
      rad = 180.0D0/pi
!   Obtain Light data
      call lidata(lam,aw,bw,ac,bc,bpic)
!   Quanta conversion
!      h = 6.6256E-34   !Plancks constant J sec
      planck = 6.6256E-34   !Plancks constant J sec
      c = 2.998E8      !speed of light m/sec
!      hc = 1.0/(h*c)
      hc = 1.0/(planck*c)
      oavo = 1.0/6.023E23   ! 1/Avogadros number
      hcoavo = hc*oavo
      do nl = 1,nlt
       rlamm = float(lam(nl))*1.0E-9  !lambda in m
       WtoQ(nl) = rlamm*hcoavo        !Watts to quanta conversion
      enddo
!   CDOM and detrital absorption exponent
      Sdom = 0.014
      Sdet = 0.013      !from Gallegos et al., 2011, small detritus nm-1
      do nl = 1,nlt
       rlam = float(lam(nl))
       excdom(nl) = exp(-Sdom*(rlam-443.0))
       exdet(nl)  = exp(-Sdet*(rlam-440.0))
      enddo

!  Spectral reflectance factor (from Frouin et al., 1996)
      do nl = 1,nlt
       rlam = float(lam(nl))
       if (lam(nl) .lt. 900)then
        t = exp(-(aw(nl)+0.5*bw(nl)))
        tlog = alog(1.0E-36+t)
        fac = a0 + a1*tlog + a2*tlog*tlog + a3*tlog*tlog*tlog
        wfac(nl) = min(fac,1.0)
        wfac(nl) = max(fac,0.0)
       else
        fac = b0 + b1*rlam + b2*rlam*rlam + b3*rlam*rlam*rlam
        wfac(nl) = max(fac,0.0)
       endif
      enddo

      return
      end
