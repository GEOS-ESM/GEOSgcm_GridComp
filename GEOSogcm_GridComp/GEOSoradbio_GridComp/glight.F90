      subroutine glight(km,IS_MIDNIGHT,cosz,lam,aw,bw,ac,bc,bpic,      &
       excdom,exdet,WtoQ,Ed,Es,H,Phyto,cdet,pic,cdc,tirrq,cdomabsq,    &
       avgq,dt)
!       deltaE)

#include "definebio.h"
#include "comlte.h"
      integer lam(nlt)
      real Ed(nlt),Es(nlt)
      real WtoQ(nlt)
      real H(km)
      real Phyto(km,nchl)
      real P(km,ntyp)
      real cdet(km),pic(km),cdc(km)
      real excdom(nlt),exdet(nlt)
      real tirrq(km),cdomabsq(km),avgq(km)
      real deltaE(km)
      real aw(nlt),bw(nlt)
      real ac(nchl,nlt),bc(nchl,nlt)
      real bpic(nlt)
      logical IS_MIDNIGHT
      data rn /1.341/  !refractive index of seawater

      P = 0.0
      P(:,nnut+1) = Phyto(:,1)
      P(:,nnut+2) = Phyto(:,2)
      P(:,nnut+3) = Phyto(:,3)
      P(:,nnut+4) = Phyto(:,4)
      P(:,nnut+5) = Phyto(:,5)
      P(:,nnut+6) = Phyto(:,6)
      P(:,nds)    = cdet
      P(:,ncs+3)  = pic
      P(:,ncs+4)  = cdc

!  Set daily parameters for ocean irradiance.
      if (IS_MIDNIGHT)then
       call daysetrad(km,avgq,dt)
      endif

!!  Compute acdom
!      do k = 1,km
!       actot450 = 0.0
!       atot450 = 0.0
!       if (H(k) < 1.0E10)then
!        do n = 1,nchl
!         actot450 = actot450 + P(k,nnut+n)*ac(n,nl450)
!!  from Morel and Maritorena, 2002
!        enddo
!        atot450 = aw(nl450) + actot450
!        do nl = 1,nlt
!         acdom(k,nl) = 0.2*atot450*excdom(nl)
!        enddo
!       endif
!      enddo

!  Compute average cosine for direct irradiance in the water 
!  column given solar zenith angle (in radians) at surface.
      sza = acos(cosz)
      sinszaw = sin(sza)/rn
      szaw = asin(sinszaw)
      rmudl = 1.0/cos(szaw)   !avg cosine direct (1 over)
      rmud = min(rmudl,1.5)
      rmud = max(rmud,0.0)

!  Compute underwater irradiance and total quanta
      call edeu(km,lam,aw,bw,ac,bc,bpic,WtoQ,Ed,Es,H,P,         &
       excdom,exdet,rmud,tirrq,cdomabsq,avgq)
!       deltaE)

      return
      end
