      subroutine setbio(rmumax,cnratio,cfratio,solFe,drate,            &
       Rm,rlamz,greff,dratez1,dratez2,regen,remin,rkn,rks,rkf,         &
       rik,cchl,wss,axs,                                               &
       Pdeep)
!  Computes constants over entire run of model, reads in required
!  data files, and otherwise obtains one-time-only information 
!  necessary for the run.
!  Biological constituents
!  P(1) = nitrate (uM)
!  P(2) = ammonium (uM)
!  P(3) = silica (uM)
!  P(4) = iron (nM)
!  P(5) = diatoms (mg chl m-3)
!  P(6) = chlorophytes (mg chl m-3)
!  P(7) = cyanobacteria (mg chl m-3)
!  P(8) = coccolithophores (mg chl m-3)
!  P(9) = dinoflagellates (mg chl m-3)
!  P(10) = phaeocystis (mg chl m-3)
!  P(11) = herbivores (mg chl m-3)
!  Detrital components
!  P(12) = N/C detritus (ugC l-1)
!  P(13) = silica detritus (uM)
!  P(14) = iron detritus (nM)
!  Carbon components
!  P(15) = DOC (uM)
!  P(16) = DIC (uM)
!  P(17) = Alkalinity  uM(C)
!  P(18) = PIC  (ugC/l)
!  P(19) = CDC  (uM)
!  pCO2 (uatm)

#include "definebio.h"
      real, parameter :: sday=86400.0    !seconds per day
      real rmumax(nchl)     !max phyto growth rate at 20oC, /d
      real rik(3,nchl)      !light saturation parameter umol quanta/m2/s
      real cchl(3)          !C:chl ratio g:g for 3 light states
!      real fescavrate(2)    !Fe scavenging rate /s
      real remin(ndet)      !detrital remineralization rate /s
      real axs(ntyp,2)      !detrital sinking a and b coefficients
      real rkn(nchl),rks(nchl) !half-sat constants for N and S (uM)
      real rkf(nchl)        !half-sat constant for Fe (nM)
      real wsd(ntyp)        !particle sinking rate m/d
      real wss(ntyp)        !particle sinking rate m/s
      real Pdeep(ntyp)

      pHsfc = 8.0
      pHmin = 7.5
      pHmax = 8.6

!  Physical parameters needed for biolog
      solFe = 0.02          !solubility of iron

!  Zooplankton/grazing parameters
      rlamz   = 1.0           !Ivlev constant
      greff   = 0.25          !grazing efficiency
      greff   = 0.10          !grazing efficiency
      dratez1 = 0.1/sday      !zooplankton death rate/s
      dratez2 = 0.5/sday      !zooplankton death rate/s
      regen   = 0.25          !regeneration fraction
      regen   = 0.10          !regeneration fraction
      Rm      = 1.20/sday     !max zoopl. growth rate/s
!      Rm      = 1.40/sday     !max zoopl. growth rate/s

!  Phytoplankton group parameters
      drate  = 0.05/sday      !phytoplankton death rate/s
      rkn    = 0.0
      rks    = 0.0
      rkf    = 0.0
      wsd    = 0.0
      wss    = 0.0
      rmumax = 0.0
      rik    = 0.0
      Pdeep    =  0.0   !bottom BC 
      Pdeep(1) = 32.0   !bottom BC for nitrate
      Pdeep(2) =  0.1   !bottom BC for ammonium
      Pdeep(3) = 60.0   !bottom BC for silica
      Pdeep(4) =  0.6   !bottom BC for iron from Archer and Johnson 2000
      Pdeep(nds:nde)  =  0.0   !bottom BC for detritus
      Pdeep(ncs)   = 0.0  !bottom BC for DOC
      Pdeep(ncs+1)  = 2345.0    !bottom BC for DIC uM from extrapolation
                                !of last layer
      Pdeep(ncs+2) = 2445.0     !bottom BC for Alk uM from extrapolation
                                !of last layer
      Pdeep(ncs+3) = 0.0        !bottom BC for PIC ug l-1

!  Carbon:chl ratios for different adaptation states
      cchl(1) = 25.0
      cchl(2) = 50.0
      cchl(3) = 80.0
!      cchl(1) = 30.0
!      cchl(2) = 60.0
!      cchl(3) = 100.0
!      cchl(1) = 20.0
!      cchl(2) = 60.0
!      cchl(3) = 100.0
      cnratio = 106.0/16.0*12.0    !C:N ratio (ugl:uM)
      csratio = 106.0/16.0*12.0    !C:Si ratio (ugl:uM)
      cfratio = 150000.0*12.0*1.0E-3    !C:Fe ratio (ugl:nM)

#if NCHL_DEFINED > 0
!  Diatoms
      n = 1
!      rmumax(n)   = 2.50     !u max in /day at 20C
      rmumax(n)   = 3.00     !u max in /day at 20C
      wsd(n+nnut) = 0.75     !sinking rate in m/day
      rik(1,n)    = 90.0     !low light-adapted Ik (<50 uE/m2/s)
      rik(1,n)    = 90.0/3.0     !low light-adapted Ik (<50 uE/m2/s)
      rik(2,n)    = 93.0     !medium light-adapted Ik (50-200 uE/m2/s)
      rik(3,n)    = 184.0    !high light adapted Ik (>200 uE/m2/s)
      rkn(n)      = 1.0      !M-M half-sat constant for nitrogen
      rks(n)      = 0.2      !M-M half-sat constant for silica
      rkf(n)      = 0.12     !M-M half-sat constant for iron
#endif
#if NCHL_DEFINED > 1
!  Chlorophytes
      n = 2
      rmumax(n)   = rmumax(n-1)*0.840
      rmumax(n)   = rmumax(n-1)*0.8
      wsd(n+nnut) = 0.25
      rik(1,n)    = rik(1,n-1)*1.077
      rik(2,n)    = rik(2,n-1)*0.935
      rik(3,n)    = rik(3,n-1)*0.781
      rkn(n)      = rkn(n-1)*0.75
!      rkn(n)      = rkn(n-1)*0.667  !1/3 distance bet. cocco and dia
      rkf(n)      = rkf(n-1)*0.835  !midway between cocco's and diatoms
!      rkf(n)      = rkf(n-1)*0.779  !1/3 distance bet. cocco and dia
#endif
#if NCHL_DEFINED > 2
!  Cyanobacteria
      n = 3
      rmumax(n)   = rmumax(n-2)*0.670
      rmumax(n)   = rmumax(n-2)*0.4
      wsd(n+nnut) = 0.0085
      rik(1,n)    = rik(1,n-2)*0.723
      rik(2,n)    = rik(2,n-2)*0.710
      rik(3,n)    = rik(3,n-2)*0.256
      rkn(n)      = rkn(n-2)*0.45
      rkn(n)      = rkn(n-2)*0.20
      rkf(n)      = rkf(n-2)*0.67  !equals cocco 
#endif
#if NCHL_DEFINED > 3
!  Coccolithophores
      n = 4
!      rmumax(n)   = rmumax(n-3)*0.663  !all coccos
!      rmumax(n)   = rmumax(n-3)*0.755  !E. huxleyi only
      rmumax(n)  = rmumax(n-3)*0.781  !E. huxleyi only (no Sunda/Hunts)
      wsd(n+nnut) = 0.82
      wsd(n+nnut) = 0.648
      rik(1,n)    = rik(1,n-3)*0.623
      rik(2,n)    = rik(2,n-3)*0.766
      rik(3,n)    = rik(3,n-3)*0.899
      rkn(n)      = rkn(n-3)*0.5
      rkf(n)      = rkf(n-3)*0.67
#endif
#if NCHL_DEFINED > 4
!  Dinoflagellates
      n = 5
      rmumax(n)   = rmumax(n-4)*0.335
      rmumax(n)   = rmumax(n-4)*0.935
      wsd(n+nnut) = 0.0
      rik(1,n)    = rik(1,n-4)*1.321
      rik(2,n)    = rik(2,n-4)*1.381
      rik(3,n)    = rik(3,n-4)*1.463
      rkn(n)      = rkn(n-4)*1.0
      rkf(n)      = rkf(n-4)*0.67
#endif
#if NCHL_DEFINED > 5
!  Phaeocystis
      n = 6
!      rmumax(n)   = rmumax(n-5)*0.655
!      rmumax(n)   = 1.3         !Schoemann et al 2005 for all Pha spp
!      rmumax(n)   = 1.0         !Schoemann et al 2005
      rmumax(n)   = 1.1
      wsd(n+nnut) = 0.75         !Peperzak et al., 2003
      rik(1,n)    = rik(1,n-5)*2.219   !Van Hilst and Smith 2002
      rik(2,n)    = rik(2,n-5)*1.263
      rik(3,n)    = rik(3,n-5)*1.064
      rkn(n)      = rkn(n-5)*0.9
      rkf(n)      = rkf(n-5)*0.1    !Arrigo et al. 1998, 2003
      rkf(n)      = rkf(n-5)*2.0    !Coale et al. 2003
#endif
!  Detrital sinking rates m/d
!      wsd(nds)   = 10.0     !nitrogen
!      wsd(nds+1) = 1.0     !silica
!      wsd(nde)   = 0.3     !iron
      axs = 0.0
      axs(nds,1)   = 1.0    !nitrogen a coef
      axs(nds,2)   = 0.1    !nitrogen b coef
      axs(nds+1,1) = 1.0    !silica a coef
      axs(nds+1,2) = 5.0    !silica b coef
      axs(nde,1)   = 0.1    !iron a coef
      axs(nde,2)   = 2.0    !iron b coef

!  PIC sinking rate m/d
      wsd(ncs+3) = 7.0
      axs(ncs+3,1) = 0.01    !PIC a coef
      axs(ncs+3,2) = 2.0    !PIC b coef

!  Diatom exponential sinking rate
      axs(nnut+1,1) = 0.01    !diatom a coef
      axs(nnut+1,2) = 5.0    !diatom b coef

      wss = wsd/sday  !convert to m/s

!  Detrital remineralization rates /s
      remin(1) = 0.01/sday             !nitrogen 
      remin(2) = 0.001/sday            !silica 
      remin(3) = 0.5/sday             !iron 
!      fescavrate(1) = 2.74E-5/sday    !low fe scavenging rate/s
!      fescavrate(2) = 50.0*fescavrate(1) !high fe scavenging rate/s

      return
      end
