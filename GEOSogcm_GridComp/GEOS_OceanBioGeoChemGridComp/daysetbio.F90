      subroutine daysetbio(km,rik,cchl,cchlratio,                   &
                           avgq,gcmax,P,H,rikd,wssc)
 
!  Sets daily parameters for bio.
 
#include "definebio.h"
      real, parameter :: sday=86400.0   !seconds per day
      real rik(3,nchl)      !light saturation parameter umol quanta/m2/s
      real cchl(3)          !C:chl ratio g:g for 3 light states
!      real fescavrate(2)         !Fe scavenging rate /s
      real T(km),avgq(km),gcmax(km),H(km)
      real P(km,ntyp)
      real rikd(km,nchl)
      real wssc(km)
      real cchlratio(km)

!  Find bottom of grid
      k = 1
      do while (k .le. km )
       if (H(k) >= 1.0E10) exit
       kend = k
       k = k+1
      enddo

!!  Compute daily variables
!      tfac20 = 0.81*exp(0.0631*20.0)    !Bissinger et al., 2008
!      do k = 1,kend
!!   Compute new max growth; use Bissinger et al,. 2008,
!!    convert to /s units (instead of /day) by dividing by
!!    86400.  Normalized to the max growth rate of diatoms
!       tfact = 0.81*exp(0.0631*T(k))    !Bissinger et al., 2008
!       tfac(k) = tfact/tfac20
!       do n = 1,nchl
!        rmuplsr(k,n) = (rmumax(n)*tfac(k))/sday 
!       enddo
!#if NCHL_DEFINED > 2
!!  Additional T-dependent factor for cyanobacteria
!       n = 3
!       tfac2 = T(k)*0.029411 + 0.55823
!       tfac2 = min(tfac2,1.0)
!       rmuplsr(k,n) = tfac2*rmuplsr(k,n)
!#endif
!!  Above provides for lower growth of cyanobacteria in cold water.
!!  The principle is to keep the delta t difference between cyano to 
!!  diatoms constant at low T
!      enddo
 
!!  Zooplankton grazing temperature dependence
!!  do not divide by 24: it is a factor
!      k = 1
!      tzoo = 0.06*exp(0.1*T(k)) + 0.70

!  Set N:chl ratio and photoadaptation
      rlightlo = 20.0  !threshold for low light adaptation (quanta/m2/s)
      rlighthi = 175.0 !threshold for high light adaptation 
      rmidq = (rlighthi-rlightlo)*0.5 + rlightlo
      do k = 1,kend
       if (avgq(k) .gt. 0.0)then
        if (avgq(k) .gt. rlighthi)then
         ikd = 3
         cchld = cchl(ikd)
         do n = 1,nchl
          rikd(k,n) = rik(ikd,n)
         enddo
        else if (avgq(k) .lt. rlightlo)then
         ikd = 1
         cchld = cchl(ikd)
         do n = 1,nchl
          rikd(k,n) = rik(ikd,n)
         enddo
        else
         ikd = 2
         if (avgq(k) .gt. rmidq)then
          fac = (avgq(k)-rmidq)/(rlighthi-rmidq)
          cchld = (1.0-fac)*cchl(2) + fac*cchl(3)
          do n = 1,nchl
           rikd(k,n) = (1.0-fac)*rik(2,n) + fac*rik(3,n)
          enddo
         else
          fac = (avgq(k)-rlightlo)/(rmidq-rlightlo)
          cchld = (1.0-fac)*cchl(1) + fac*cchl(2)
          do n = 1,nchl
           rikd(k,n) = (1.0-fac)*rik(1,n) + fac*rik(2,n)
          enddo
         endif
        endif
       else
        ikd = 1
        cchld = cchl(ikd)
        do n = 1,nchl
         rikd(k,n) = rik(ikd,n)
        enddo
       endif
       do n = nnut+1,npe+nzoo
        P(k,n) = P(k,n)*cchlratio(k)/cchld
       enddo
       cchlratio(k) = cchld
      enddo
      avgq = 0.0
 
#if NCHL_DEFINED > 3
!  Adjustable sinking rate for cocco's: range = 0.3 to 1.4 m/day
      do k = 1,kend
       gcmaxd = gcmax(k)*sday
       wstmp = 0.752*gcmaxd + 0.225
!       wstmp = 0.1*exp(1.4*gcmaxd)
       wsdc = max(wstmp,0.01)
       wsdc = min(wsdc,1.4)
       wssc(k) = wsdc/sday
      enddo
      do k = 1,kend
       gcmax(k) = 0.0
      enddo
#endif
!!  Fe scavenging 
!      do k = 1,kend
!       if (P(k,4) .lt. 0.6)then
!        fescav(k) = fescavrate(1)*P(k,4)
!       else
!        fescav(k) = fescavrate(1)*P(k,4)+fescavrate(2)*(P(k,4)-0.6)
!       endif
!      enddo
 
      return
      end
