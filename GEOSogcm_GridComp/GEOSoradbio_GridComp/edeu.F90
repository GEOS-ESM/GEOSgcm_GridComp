      subroutine edeu(km,lam,aw,bw,ac,bc,bpic,WtoQ,Ed,Es,H,P,    &
       excdom,exdet,rmud,tirrq,cdomabsq,avgq)

!  Model of irradiance in the water column.  Accounts for three 
!  irradiance streams:

!  Edz = direct downwelling irradiance
!  Esz = diffuse downwelling irradiance
!  Euz = diffuse upwelling irradiance

!  Propagation is done in energy units, tests are done in quanta,
!  final is quanta for phytoplankton growth.

#include "definebio.h"
#include "comlte.h"
      integer lam(nlt)
      real Ed(nlt),Es(nlt)
      real WtoQ(nlt)
      real H(km)
      real P(km,ntyp)
      real acdom(nlt)
      real excdom(nlt),exdet(nlt)
      real tirrq(km),cdomabsq(km),avgq(km)
      real deltaE(km)
      real aw(nlt),bw(nlt)
      real ac(nchl,nlt),bc(nchl,nlt)
      real bpic(nlt)
      real Edz(nlt,km),Esz(nlt,km)
      real Euz(nlt,km)
      real Edtop(nlt),Estop(nlt)
      real fchl(nchl)
      real bbrc(nchl)
      data bbrc /0.002, 0.00071, 0.0032, 0.00071, 0.0029, 0.002/
      data bbrw /0.5/         !backscattering to total scattering ratio
      data bbrpic /0.01/      !Balch et al., 1996, low end
      data bbrd /0.005/       !Gallegos et al., 2011 small detritus
      data adstar /8.0E-5/    !Gallegos et al., 2011 m2/mg
      data bdstar /0.00115/   !Gallegos et al., 2011 m2/mg
      data acdomstar /2.98E-4/  !Yacobi et al., 2003 m2/mg
      data Dmax /500.0/       !depth at which Ed = 0
      data npst,npnd /3,17/

!  Constants and initialize
      rmus = 1.0/0.83            !avg cosine diffuse down
      tirrq = 0.0
      cdomabsq = 0.0
      deltaE = 0.0
      acdom = 0.0

      Ebot = 0.0
      do nl = 1,nlt
!      do nl = npst,npnd
       Edtop(nl) = Ed(nl)
       Estop(nl) = Es(nl)
       Ebot = Ebot + (Ed(nl)+Es(nl))
      enddo
!  Convert to quanta: divide by Avos # to get moles quanta; then mult by
!  1E6 to get uM or uEin
      Ebotq = 0.0
      do nl = 1,nlt
!      do nl = npst,npnd   !PAR range only 350-700nm
       Ebotq = Ebotq + (Edtop(nl)+Estop(nl))*WtoQ(nl)*1.0E6
      enddo
      do k = 1,km
       if (H(k) < 1.0E10)then
        Etop = Ebot
        Etopq = Ebotq
        zd = min(Dmax,H(k))
        zirr = 0.0
        zirrq = 0.0
        do nl = 1,nlt
!         do nl = npst,npnd
         Edz(nl,k) = 0.0
         Esz(nl,k) = 0.0
         Euz(nl,k) = 0.0
         actot = 0.0
         bctot = 0.0
         bbctot = 0.0
         bbctot = 0.0
         Plte = max(P(k,ncs+4),0.0)
         acdom(nl) = Plte*12.0*acdomstar*excdom(nl) !cdoc in units of uM
!   No detritus or PIC optics
         do n = 1,nchl
          Plte = max(P(k,nnut+n),0.0)
          actot  = actot  + Plte*ac(n,nl)
          bctot  = bctot  + Plte*bc(n,nl)
          bbctot = bbctot + Plte*bbrc(n)*bc(n,nl)
         enddo
!         a  = aw(nl) + acdom(nl) + actot
!         bt = bw(nl) + bctot
!         bb = bbrw*bw(nl) + bbctot
!   Detritus optics
         Plte = max(P(k,nds),0.0)
         adet = Plte*adstar*exdet(nl)
         bdet = Plte*bdstar*(555.0/float(lam(nl))**0.5)
         a  = aw(nl) + acdom(nl) + actot + adet
!         bt = bw(nl) + bctot + bdet
!         bb = bbrw*bw(nl) + bbctot + bbrd*bdet
!   Detritus and PIC scattering
         Plte3 = max(P(k,ncs+3),0.0)
         bt = bw(nl) + bctot + bdet + bpic(nl)*Plte3
         bb = bbrw*bw(nl) + bbctot + bbrd*bdet + bbrpic*bpic(nl)*Plte3
         bb = max(bb,0.0002)
         if (Edtop(nl) .ge. 1.0E-4 .or. Estop(nl) .ge. 1.0E-4)then
          call radmod(zd,Edtop(nl),Estop(nl),rmud,a,bt,bb,           &
            Dmax,Edz(nl,k),Esz(nl,k),Euz(nl,k))
         endif
         Edtop(nl) = Edz(nl,k)
         Estop(nl) = Esz(nl,k)
         zirr = zirr + (Edz(nl,k)+Esz(nl,k)+Euz(nl,k))
        enddo
        Ebot = zirr
        deltaE(k) = Etop - Ebot
        do nl = 1,nlt
!        do nl = npst,npnd
         sumQ = (Edz(nl,k)+Esz(nl,k)+Euz(nl,k))*WtoQ(nl)*1.0E6
         cdomabsq(k) = cdomabsq(k) + acdom(nl)*sumQ
         zirrq = zirrq + sumQ
        enddo
        Ebotq = zirrq

!  tirrq is the average of the natural log between Etopq and Ebotq
!  converted to sqrt using the identity relationships for speedup

        tirrq(k) = sqrt(Etopq*Ebotq)*rmus
       endif
      enddo

!  Irradiance summary loops
      do k = 1,km
       avgq(k) = avgq(k) + tirrq(k)
      enddo

      return
      end
