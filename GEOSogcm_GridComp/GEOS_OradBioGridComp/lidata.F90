      subroutine lidata(lam,aw,bw,ac,bc,bpic)
 
!  Reads in radiative transfer data: specifically 
!  water data (seawater absorption and total scattering coefficients,
!  and chl-specific absorption and total scattering data for 
!  several phytoplankton groups).  PAR (350-700) begins at index 3, 
!  and ends at index 17.
 
#include "definebio.h"
#include "comlte.h"
      character*80 title
      character*80 cfle
      character cacbc*12,cabw*15,cpic*13
      real saw,sbw,sac,sbc
      character*4 cdir
      integer lam(nlt)
      real aw(nlt),bw(nlt)
      real ac(nchl,nlt),bc(nchl,nlt)
      real bpic(nlt)
      data cdir /'bcs/'/
      data cacbc,cabw,cpic /'acbc25_6.dat','abw25_morel.dat',          &
                            'pic_sigma.dat'/
 
!  Water data files
      cfle = cdir//cabw
      cfle = cabw
      open(4,file=cfle,status='old',form='formatted')
      do i = 1,6
       read(4,'(a80)')title
!       write(6,'(a80)')title
      enddo
      do nl = 1,nlt
       read(4,20)lambda,saw,sbw
!       write(6,20)lambda,saw,sbw
       lam(nl) = lambda
       aw(nl) = saw
       bw(nl) = sbw
      enddo
      close(4)
20    format(i5,f15.4,f10.4)
!  Modify for Lee et al. 2015
      aw(3)  = 0.0071
      aw(4)  = 0.0041
      aw(5)  = 0.0034
      aw(6)  = 0.0033
      aw(7)  = 0.0068
      aw(8)  = 0.0099
      aw(9)  = 0.0187
      aw(10) = 0.0400
      aw(11) = 0.0652
 
!  Phytoplankton group chl-specific absorption and total scattering 
!  data.  Chl-specific absorption data is normalized to 440 nm; convert
!  here to actual ac*(440)
      cfle = cdir//cacbc
      cfle = cacbc
      open(4,file=cfle,status='old',form='formatted')
      do i = 1,6
       read(4,'(a80)')title
      enddo
      do n = 1,nchl
       read(4,'(a80)')title
       do nl = 1,19
        read(4,30)lambda,sac,sbc
        ac(n,nl) = sac
        bc(n,nl) = sbc
       enddo
       do nl = 20,nlt
        ac(n,nl) = 0.0
        bc(n,nl) = 0.0
       enddo
      enddo
      close(4)
30    format(i4,2f10.4)

!  PIC scattering cross section m2/ugCaCO3 (PIC)  (=PIC-specific
!  scattering coefficient (from Gordon et al., 2009)
      cfle = cdir//cpic
      cfle = cpic
      open(4,file=cfle,status='old',form='formatted')
      do i = 1,5
       read(4,'(a80)')title
      enddo
      do nl = 1,nlt
       read(4,40)lamm,bpic(nl)
      enddo
      close(4)
40    format(i4,f12.6)
 
      return
      end
