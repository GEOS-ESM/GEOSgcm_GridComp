      subroutine carbon(k,cnratio,cchlratio,bn,bf,remin,P,tfac,      &
       fnoice,tzoo,tirrq,cdomabsq,gro,picdis,                        &
       atmco2,wspd,slp,T,S,H,P_tend,pco2,cflx)

!  Computes carbon cycling.  Uses Aumont et al (2002; JGR) for
!  semi-labile DOC (because of basic similarities in model
!  formulation, specifically, presense of grazers and particulate
!  detritus).
!  ncar(1) = DOC semi-labile (uM(C))
!  ncar(2) = DIC  (uM(C))
!  ncar(3) = Alkalinity  uM(C)
!  ncar(4) = PIC   ugC/l
!  ncar(5) = CDOC  uM
 
#include "definebio.h"
      real, parameter :: sday=86400.0   !seconds per day
      real, parameter :: excp=0.05      !excretion of DOC by phyto %growth
      real, parameter :: resp=0.05      !respiration of DIC by phyto %growth
      real, parameter :: excz=0.05/sday !excretion of DOC by zoopl/s
      real, parameter :: resz=0.05/sday !respiration of DIC by zoopl/s
!      real, parameter :: rlamdoc=0.034/sday !N-dependent DOC remin/s
!      real, parameter :: rlamdoc=0.017/sday !N-dependent DOC remin/s
!      real, parameter :: rlamdoc=0.005/sday !N-dependent DOC remin/s
      real, parameter :: rlamdoc=0.010/sday !N-dependent DOC remin/s
      real, parameter :: rkdoc1=0.3*10.0 !N-dep half sat uM(PO4) modified
!                                for nitrate by mult*10.0, based on 
!                                Conkright et al. 1994
      real, parameter :: rkdoc2=15.0    !DOC-dep half-sat uM(C)
      real, parameter :: rlampoc=0.05/sday  !detrital breakdown/s
      real, parameter :: uMtomgm3=12.0  !conversion uM to mg/m3
      real, parameter :: Pzo=1.0*uMtomgm3/50.0 !zoopl half-sat for 
                                  !DOC excretion mg/m3(chl,assuming
                                  !C:chl ratio of 50))
      real, parameter :: frdoc=0.7      !fraction DOC:CDOC
      real, parameter :: awan=0.337/(3.6E+5) !piston vel coeff., from 
!                                     Wanninkhof 1992, but adjusted 
!                                     by OCMIP, and converted from 
!                                     cm/hr to m/s
      real, parameter :: stdslp=1013.25 !standard sea level pressure mb
!      real, parameter :: atmco2=371.3  !uatm or ppmv (equivalent); 
!                             global mean 2000-2003 from OCMIP
      real mgchltouMC
      real remin(ndet)      !detrital remineralization rate /s
      real docexcp(nchl),cdocexcp(nchl),dicresp(nchl)
      real P(ntyp),gro(nchl)
      real P_tend(ntyp)
      real tfac
      real tirrq,cdomabsq
      real fnoice,tzoo

      bs = 2.0*bn
      frcdoc = 1.0-frdoc
      mgchltouMC = cchlratio/uMtomgm3
!  Detrital, bacterial, and grazing components
!   DOC
      Pnpe1 = max(P(npe+1),0.0)
      rmmzoo = Pnpe1/(Pzo+Pnpe1)
!      docexcz  = frdoc*excz*rmmzoo*Pnpe1        !zoopl production DOC
!      cdocexcz = frcdoc*excz*rmmzoo*Pnpe1       !zoopl production CDOC
      docexcz  = excz*rmmzoo*Pnpe1        !zoopl production DOC
      cdocexcz = 0.0       !zoopl production CDOC
      P_tend(npe+1) = P_tend(npe+1) - (docexcz+cdocexcz)*fnoice
      P_tend(1) = P_tend(1) + bn*(docexcz+cdocexcz)*fnoice
      P_tend(4) = P_tend(4) + bf*(docexcz+cdocexcz)*fnoice
      P1 = max(P(1),0.0)
      rndep = rlamdoc*(P1/(rkdoc1+P1))
      Pncs = max(P(ncs),0.0)
      Pncs4 = max(P(ncs+4),0.0)
      Pnds = max(P(nds),0.0)
      docdep = Pncs/(rkdoc2+Pncs)
!      docbac = tfac*frdoc*rndep*docdep*Pncs   !bacterial loss DOC
!      docdet = tfac*frdoc*rlampoc*Pnds        !detrital production DOC
      docbac = tfac*rndep*docdep*Pncs   !bacterial loss DOC
      docdet = tfac*rlampoc*Pnds        !detrital production DOC
      doctend = docexcz*mgchltouMC + docdet/uMtomgm3-docbac
      P_tend(ncs) = P_tend(ncs) + doctend*fnoice
!      cdocbac = tfac*frcdoc*rndep*docdep*Pncs4  !bacterial loss CDOC
!      cdocdet = tfac*frcdoc*rlampoc*Pnds4       !detrital production CDOC
      cdocbac = 0.0  !bacterial loss CDOC
      cdocdet = 0.0       !detrital production CDOC
      cdoctend = cdocexcz*mgchltouMC + cdocdet/uMtomgm3-cdocbac
      P_tend(ncs+4) = P_tend(ncs+4) + cdoctend*fnoice
!   adjust detritus
      P_tend(nds) = P_tend(nds) - (docdet+cdocdet)*fnoice  !carbon/nitrogen detritus
!   equivalent amount of DOC from detritus goes into nitrate, bypassing
!     DON
      P_tend(1) = P_tend(1) + (docdet+cdocdet)/cnratio*fnoice
!   DIC
      dicresz = tzoo*resz*Pnpe1 !zoopl production DIC (resp)
      P_tend(npe+1) = P_tend(npe+1) - dicresz*fnoice
      P_tend(1) = P_tend(1) + bn*dicresz*fnoice
      P_tend(4) = P_tend(4) + bf*dicresz*fnoice
      P_tend(ncs+1) = (dicresz*mgchltouMC + docbac + cdocbac            &
                        + tfac*remin(1)*Pnds/uMtomgm3)*fnoice
!  Phytoplankton components related to growth
      picgro = 0.0
      if (tirrq .gt. 0.0)then
       sumdoc = 0.0
       sumutk = 0.0
       sumres = 0.0
       docexcp = 0.0
       cdocexcp = 0.0
       dicresp = 0.0
       do n = 1,nchl
        totgro = gro(n)
        docexcp(n) = frdoc*excp*totgro    !phyto production DOC
        cdocexcp(n) = frcdoc*excp*totgro  !phyto production CDOC
        dicresp(n) = resp*totgro          !phyto production DIC 
!   dont account for ice here -- it is incorporated in gro/totgro
        P_tend(n+nnut) = P_tend(n+nnut) - (docexcp(n)+dicresp(n))
        P_tend(n+nnut) = P_tend(n+nnut) - cdocexcp(n)
        P_tend(1) = P_tend(1) + bn*(docexcp(n)+cdocexcp(n)+dicresp(n))
        P_tend(4) = P_tend(4) + bf*(docexcp(n)+cdocexcp(n)+dicresp(n))
        sumdoc = sumdoc + docexcp(n)
        sumcdoc = sumcdoc + cdocexcp(n)
        sumutk = sumutk + totgro
        sumres = sumres + dicresp(n)
       enddo
       P_tend(3) = P_tend(3) + bs*(docexcp(1)+cdocexcp(1)+dicresp(1))
       P_tend(ncs) = P_tend(ncs) + sumdoc*mgchltouMC
       P_tend(ncs+4) = P_tend(ncs+4) + sumcdoc*mgchltouMC
       P_tend(ncs+1) = P_tend(ncs+1) + ((sumres-sumutk)*mgchltouMC)
       picgro = 0.25*(gro(4)-resp*gro(4)) !PIC production and respiration
      endif
!   PIC production/respiration and losses
      Pncs3 = max(P(ncs+3),0.0)
      picd = picdis*fnoice
      P_tend(ncs+3) = picgro*cchlratio - picd*Pncs3  !PIC ugC/l
      P_tend(ncs+1) = P_tend(ncs+1) - picgro*mgchltouMC & !loss DIC uM
                      + picd*Pncs3/12.0      !PIC dissolve uM
      P_tend(ncs+2) = -2.0*picgro*mgchltouMC &  !Alkalinity from PIC prod
                      +2.0*picd*Pncs3/12.0 & !from PIC dissolution
                      -P_tend(1) + P_tend(2)
!   Destruction of CDOC
      destroycdoc = 1.0E-6*cdomabsq
      P_tend(ncs+4) = P_tend(ncs+4) - destroycdoc*fnoice
      P_tend(ncs+1) = P_tend(ncs+1) + destroycdoc*fnoice  !reverts to DIC


!  Surface fluxes of carbon 
      if (k .eq. 1)then
!  pCO2
       call ppco2(slp,atmco2,T,S,P,ff,pco2)
!  Update DIC for sea-air flux of CO2
       Ts = T
       scco2 = 2073.1 - 125.62*Ts + 3.6276*Ts**2 - 0.043219*Ts**3
       scco2arg = (scco2/660.0)**(-0.5)
       wssq = wspd*wspd
!       rkwco2 = awan*(wssq+wspdvar)*scco2arg !units of m/s
       rkwco2 = awan*wssq*scco2arg !units of m/s
       ffuatm = ff*1.0E-6                    !mol/kg/uatm
       xco2 = atmco2*slp/stdslp              !uatm
       deltco2 = (xco2-pco2)*ffuatm*1024.5   !mol/m3
       flxmolm2 = rkwco2*deltco2             !units of mol/m2/s
       cflx = -flxmolm2*12.0E6*fnoice        !units of ugC/m2/s
       flxmolm3 = flxmolm2/H                 !mol/m3/s
       P_tend(ncs+1) = P_tend(ncs+1) + flxmolm3*1000.0*fnoice
                                             !uM/s
      endif

      return
      end

