      subroutine lidatatm(Fobar,thray,oza,awv,ao,aco2)
 
!  Reads in radiative transfer data: specifically 
!  water data (seawater absorption and total scattering coefficients,
!  and chl-specific absorption and total scattering data for 
!  several phytoplankton groups).  PAR (350-700) begins at index 3, 
!  and ends at index 17.
 
#include "definebio.h"
#include "comlte.h"
      character*80 title
      character*80 cfle
      character catmo*11,cacbc*11,cabw*15,cpic*13
      real saw,sbw,sac,sbc
      character*4 cdir
      real Fobar(nlt),thray(nlt),oza(nlt),awv(nlt),ao(nlt),aco2(nlt)
      data cdir /'bcs/'/
      data catmo /'atmo25b.dat'/

!  Atmospheric data file
      cfle = cdir//catmo
      cfle = catmo
      open(4,file=cfle,status='old',form='formatted')
      read(4,'(a80)')title
      read(4,'(a80)')title
      read(4,'(a80)')title
      read(4,'(a80)')title
      do nl = 1,nlt
       read(4,10)ilam,sFobar,sthray,soza,sawv,sao,saco2
!       write(6,10)ilam,sFobar,sthray,soza,sawv,sao,saco2
       Fobar(nl) = sFobar
       thray(nl) = sthray
       oza(nl) = soza
       awv(nl) = sawv
       ao(nl) = sao
       aco2(nl) = saco2
      enddo
      close(4)

10    format(i5,6f11.4)
      return
      end
