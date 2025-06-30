      subroutine sink(km,kend,DH,ws,P,cnratio,P_tend)
 
!  Sinking of biological particles.  Upwind differencing.
 
#include "definebio.h"
      real, parameter :: uMtomgm3=12.0  !conversion uM to mg/m3
      real DH(km)
      real P_tend(km,ntyp)
      real P(km,ntyp)
      real ws(km,ntyp)
 
      area = 1.0            !for convenience area = 1.0m2
 
!  Convert all variables to mg/s, convert ws to m3/s
      do k = 1,kend
       do n = 1,ntyp
        P_tend(k,n) = P_tend(k,n)*DH(k)*area
        ws(k,n) = ws(k,n)*area
       enddo
      enddo
 
!  Sinking as mg/s
      do k = 1,kend-1
       do n = nnut+1,npe
        trnd = P(k,n)*ws(k,n)
        P_tend(k,n)   = P_tend(k,n)   - trnd
        P_tend(k+1,n) = P_tend(k+1,n) + trnd
       enddo  
       do n = nds,nde
        trnd = P(k,n)*ws(k,n)
        P_tend(k,n)   = P_tend(k,n)   - trnd
        P_tend(k+1,n) = P_tend(k+1,n) + trnd
       enddo  
       n = ncs+3
       trnd = P(k,n)*ws(k,n)
       P_tend(k,n)   = P_tend(k,n)   - trnd
       P_tend(k+1,n) = P_tend(k+1,n) + trnd
      enddo 

!  Let detritus in the last level sink to the bottom and be immediately
!  remineralized to nutrients.  
      k = kend
!   nitrate detritus and nitrate (also carbon)
      n = nds
      trnd = P(k,n)*ws(k,n)
      P_tend(k,n) = P_tend(k,n) - trnd
      P_tend(k,1) = P_tend(k,1) + trnd/cnratio
      P_tend(k,ncs) = P_tend(k,ncs) + trnd/uMtomgm3
!   silica and iron detritus
      do n = nds+1,nde
       trnd = P(k,n)*ws(k,n)
       P_tend(k,n)   = P_tend(k,n) - trnd
       P_tend(k,n-10) = P_tend(k,n-10) + trnd
      enddo
!  PIC is lost to sediments
      n = ncs+3
      trnd = P(k,n)*ws(k,n)
      P_tend(k,n) = P_tend(k,n) - trnd

!  Convert to concentrations per second, e.g., mg/m3/s
      do k = 1,kend
       do n = 1,ntyp
        P_tend(k,n) = P_tend(k,n)/(area*DH(k))
       enddo
      enddo
 
      return
      end

