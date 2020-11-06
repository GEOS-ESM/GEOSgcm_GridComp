      subroutine setsfclte(rad,pi2,lam,Fobar,thray,oza,awv,ao,aco2,   &
       am,Vi,asl,bsl,csl,dsl,esl,fsl,ica)

!  Computes constants for global irradiance calculations, reads in 
!  required data files, and otherwise obtains one-time-only information 
!  necessary for the run.

#include "definebio.h"
#include "comlte.h"
      integer lam(nlt)
      real Fobar(nlt),thray(nlt),oza(nlt),awv(nlt),ao(nlt),aco2(nlt)
      integer ica(nlt)
      real asl(ncld),bsl(ncld),csl(ncld),dsl(ncld),esl(ncld),fsl(ncld)

!  Degrees to radians conversion
      pi = dacos(-1.0D0)
      pi2 = pi*2.0
      rad = 180.0D0/pi

!  Obtain Light data
      call lidatatm(Fobar,thray,oza,awv,ao,aco2)
      am = 1.0
      Vi = 25.0
      call rdslingo(lam,asl,bsl,csl,dsl,esl,fsl,ica)

      return
      end

! **********************************************************************
      subroutine rdslingo(lam,asl,bsl,csl,dsl,esl,fsl,ica)

!  Reads cloud parameters by Slingo (1989).

#include "comlte.h"
      character*50 title
      integer lam(nlt)
      integer ica(nlt)
      real asl(ncld),bsl(ncld),csl(ncld),dsl(ncld),esl(ncld),fsl(ncld)
      real rnl1(ncld),rnl2(ncld)
      real*4 rn1,rn2,a4,b4,c4,d4,e4,f4

      open(15,file='slingo.dat',status='old',form='formatted')
      do n = 1,3
       read(15,'(a50)')title
      enddo
      do n = 1,ncld
       read(15,10)rn1,rn2,a4,b4,e4,f4,c4,d4
       rnl1(n) = rn1
       rnl2(n) = rn2
       asl(n) = a4*0.01
       bsl(n) = b4
       csl(n) = c4
       dsl(n) = d4
       esl(n) = e4
       fsl(n) = f4*0.001
      enddo
      close(15)

!   Indices to relate cloud parameters to other light parameters
       do nl = 1,nlt
        do nc = 1,ncld
         lamcld = nint(rnl2(nc)*1000.0)
         if (lam(nl) .lt. lamcld)then
          ica(nl) = nc
          go to 5
         endif
        enddo
5       continue
       enddo

10    format(2f5.2,3x,2f6.3,2f6.3,1pe9.2,1pe8.2)
      return
      end
