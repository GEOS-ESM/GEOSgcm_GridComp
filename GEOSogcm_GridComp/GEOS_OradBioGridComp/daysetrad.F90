      subroutine daysetrad(km,avgq,dt)
!  Sets daily parameters for ocean irradiance.
 
#include "definebio.h"
#include "comlte.h"
      real avgq(km)
      real, parameter :: sday=86400.0  !seconds per day
 
!  Compute average quanta; Initialize light history arrays
      do k = 1,km
       avgq(k) = avgq(k)*sday/dt
      enddo
 
      return
      end

