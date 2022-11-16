      subroutine ptend(k,dt,solFe,remin,rkn,rks,rkf,cnratio,cfratio,   &
       cchlratio,wss,wssc,drate,Rm,rlamz,greff,dratez1,dratez2,        &
       regen,axs,                                                      &
       P,H,atmFe,fnoice,rikd,rmumax,gcmax,                             &
       tirrq,cdomabsq,atmco2,wspd,slp,T,S,                             &
       P_tend,gronfix,ws,Hpst,picdis,pco2,cflx,pp2)

!  Computes tendencies of biological particles.  Outputs as 
!  units/h.
!  P(1)  = nitrate (uM)
!  P(2)  = ammonium (uM)
!  P(3)  = silica (uM)
!  P(4)  = iron (nM)
!  P(5)  = diatoms (mg chl m-3)
!  P(6)  = chlorophytes (mg chl m-3)
!  P(7)  = cyanobacteria (mg chl m-3)
!  P(8)  = coccolithophores (mg chl m-3)
!  P(9)  = dinoflagellates (mg chl m-3)
!  P(10) = phaeocystis (mg chl m-3)
!  P(11) = herbivores (mg chl m-3)
!  P(12) = nitrate/carbon detritus (ugC l-1)
!  P(13) = silica detritus (uM)
!  P(14) = iron detritus (nM)
!  P(15) = DOC (uM)
!  P(16) = DIC (uM)
!  P(17) = Alkalinity  uM(C)
!  P(18) = PIC  (ugC/l)
!  P(19) = CDC  (uM)

#include "definebio.h"
      real remin(ndet)      !detrital remineralization rate /s
      real rkn(nchl),rks(nchl) !half-sat constants for N and S (uM)
      real rkf(nchl)        !half-sat constant for Fe (nM)
      real dphy(ntyp)       !phytoplankton death rate * P
      real zoo(ntyp)        !ingestion fraction by phyto group
      real rmu4(nchl),rmu3(nchl) !uptake of nutrients
      real wss(ntyp)        !particle sinking rate m/s
      real axs(ntyp,2)      !detrital sinking a and b coefs
      real P(ntyp)
      real H,tfac,fescav,wssc
      real atmFe,fnoice,tzoo
      real rikd(nchl)
      real rmuplsr(nchl)    !realized growth+respiration rate /s
      real rmumax(nchl)     !max phyto growth rate at 20oC, /d
      real gcmax,tirrq,cdomabsq
      real T,S
      real P_tend(ntyp)
      real ws(ntyp)
      real Hpst
      real pp2(nchl)
      real gro(nchl)
      real gronfix
      real viscfac
      real wsday
      real, parameter :: sday=86400.0   !seconds per day
      real, parameter :: excp=0.05   !excretion of DOC by phyto %growth
      real, parameter :: resp=0.05   !respiration of DIC by phyto %growth
      real, parameter :: frdoc=0.5   !fraction of DOC:CDOC
      real, parameter :: phygross=1.0-(excp+resp) !factor to derive net PP
      real, parameter :: uMtomgm3=12.0  !conversion uM to mg/m3
!      real, parameter :: phaeoTopt=16.3  !Phaeocystis Topt oC
      real, parameter :: phaeoTopt=1.0  !Phaeocystis Topt oC
      real, parameter :: phaeodT=13.7    !Phaeocystis dT oC
      real mgchltouMC

!  Initialize
!      P_tend = 0.0
      gro     = 0.0
      gronfix = 0.0
      ws      = 0.0
      pp2     = 0.0
      bn      = cchlratio/cnratio
      bs      = 2.0*bn
      bf      = cchlratio/cfratio
      mgchltouMC = cchlratio/uMtomgm3
      frcdoc = 1.0-frdoc

!  Compute temperature dependences
      tfac20 = 0.81*exp(0.0631*20.0)    !Bissinger et al., 2008
!   Compute new max growth; use Bissinger et al,. 2008,
!    convert to /s units (instead of /day) by dividing by
!    86400.  Normalized to the max growth rate of diatoms
       tfact = 0.81*exp(0.0631*T)    !Bissinger et al., 2008
       tfac = tfact/tfac20
       do n = 1,nchl
        rmuplsr(n) = (rmumax(n)*tfac)/sday
       enddo
!  Additional T-dependent factor for cyanobacteria
!  Provides for lower growth of cyanobacteria in cold water.
!  The principle is to keep the delta t difference between cyano to
!  diatoms constant at low T
       n = 3
       tfac2 = T*0.029411 + 0.55823
       tfac2 = min(tfac2,1.0)
       rmuplsr(n) = tfac2*rmuplsr(n)
!  Modify chloro max growth for lower growth at low T and higher at
!  high T
       n = 2
       tfac20c = 0.4*exp(0.075*20.0)
       tfactc = 0.4*exp(0.075*T)
       tfac = tfactc/tfac20c
       rmuplsr(n) = (rmumax(n)*tfac)/sday
#if NCHL_DEFINED > 5
!  Phaeocystis growth has Gaussian sensitivity to T (Schoemann et al., 
!  2005)
       n = 6
       pexarg = ((T-phaeoTopt)*(T-phaeoTopt))/(phaeodT*phaeodT)
       rmuplsrday = rmumax(n)*exp(-pexarg)
       rmuplsr(n) = rmuplsrday/sday
#endif

!  Zooplankton grazing temperature dependence
!  do not divide by 24: it is a factor
      tzoo = 0.06*exp(0.1*T) + 0.70

!  Start Model Space Loop
!   Iron + atm iron: disperse in layer and convert to nM/s
      if (k .eq. 1)then
       P_tend(4) = (atmFe*solFe)/H*0.001
      endif
!  Fe scavenging
      P4 = max(P(4),0.0)
      if (P4 .lt. 0.6)then
       fescav = (2.74E-5/sday)*P4  !low fe scavenging rate/s
      else
                !high fe scavenging rate/s
       fescav = (2.74E-5/sday)*P4+(10.0*(2.74E-5/sday))*(P4-0.3)
      endif
!   Grazing/regeneration of ammonium
      ptot = 0.0
      do n = nnut+1,npe
       ptot = ptot + max(P(n),1.0E-36)
      enddo
      ptot = max(ptot,1.0E-36)
      Pzoo = max(P(npe+1),0.0)
      gzoo = tzoo*Rm*(1.0-exp(-rlamz*ptot))*Pzoo
      dzoo1 = dratez1*Pzoo
      dzoo2 = dratez2*Pzoo*Pzoo
      P_tend(npe+1) = ((1.0-greff)*gzoo-dzoo1-dzoo2)*fnoice
      do n = nnut+1,npe
!      fraction of grazing for this group
       Pn = max(P(n),0.0)
       zoo(n) = gzoo*Pn/ptot
       dphy(n) = drate*Pn
       P_tend(n) = (-zoo(n)-dphy(n))*fnoice
      enddo
      exc = greff*gzoo
      Pnds  = max(P(nds),0.0)
      Pnds1 = max(P(nds+1),0.0)
      Pnde  = max(P(nde),0.0)
      P_tend(1) = P_tend(1) + tfac*remin(1)*Pnds/cnratio*fnoice
      P_tend(2) = bn*exc*fnoice + bn*regen*dzoo2*fnoice
      P_tend(3) = tfac*remin(2)*Pnds1*fnoice
      P_tend(4) = P_tend(4)                                          &
                        + bf*exc*fnoice                              &
                        + bf*regen*dzoo2*fnoice                      &
                        + tfac*remin(3)*Pnde*fnoice                  &
                        - fescav
      doctend = exc*mgchltouMC + regen*dzoo2*mgchltouMC
!      P_tend(ncs) = frdoc*doctend*fnoice                         !DOC
      P_tend(ncs) = doctend*fnoice                         !DOC
!      P_tend(ncs+4) = (P_tend(ncs+4)+(frcdoc*doctend))*fnoice    !CDOC
      dphyt = 0.0
      do n = nnut+1,npe
       dphyt = dphyt + dphy(n)
      enddo
!  1st detrital fraction is carbon
      P_tend(nds) = (dphyt+dzoo1)*cchlratio*fnoice                  &
                        + (1.0-regen)*dzoo2*cchlratio*fnoice        &
                        - tfac*remin(1)*Pnds*fnoice
!  2nd detrital fraction is silica
      P_tend(nds+1) = bs*dphy(nnut+1)*fnoice                        &
                        + bs*zoo(nnut+1)*fnoice                     &
                        - tfac*remin(2)*Pnds1*fnoice
!  3rd detrital fraction is iron
      P_tend(nde) = bf*(dphyt+dzoo1)*fnoice                         &
                        + bf*(1.0-regen)*dzoo2*fnoice               &
                        - tfac*remin(3)*Pnde*fnoice                 &
                        + 0.1*fescav

!   Day: Grow
      if (tirrq .gt. 0.0)then
       P1 = max(P(1),1.0E-36)     !used as a denominator
       P2 = max(P(2),1.0E-36)
       P3 = max(P(3),1.0E-36)
       P4 = max(P(4),1.0E-36)
       tirrqice = tirrq*0.01  !reduce light in ice by factor of 10
!   Light-regulated growth
#if NCHL_DEFINED > 0
!   Diatoms
       n = 1
!    Nutrient-regulated growth; Michaelis-Menton uptake kinetics
       rnut2 = P2/(rkn(n)+P2)         !ammonium
       tnit = P1/(rkn(n)+P1)          !nitrate
       tmp = 1.0 - rnut2
!    Enforce preferential utilization of ammonium
       rnut1 = min(tmp,tnit)
       rmmn = rnut1 + rnut2
       framm = rnut2/rmmn
       rmml = tirrq/(tirrq+0.5*rikd(n))
       rmmlice = tirrqice/(tirrqice+0.5*rikd(n))
!    Silica
       rmms = P3/(rks(n)+P3)
       rmmf = P4/(rkf(n)+P4)          !iron
       rlim = min(rmml,rmmn,rmms,rmmf)
       rlimice = min(rmmlice,rmmn,rmms,rmmf)
       grate = rmuplsr(n)*rlim*fnoice                                 &
                + rmuplsr(n)*rlimice*(1.0-fnoice)
       rmu4(n) = grate*framm
       rmu3(n) = grate*(1.0-framm)
       gro(n) = grate*max(P(n+nnut),0.0)
       P_tend(n+nnut) = P_tend(n+nnut) + gro(n)
!    Net primary production
       pp2(n) = gro(n)*phygross*H*cchlratio*sday
#endif
#if NCHL_DEFINED > 1
!   Chlorophytes
       n = 2
!    Nutrient-regulated growth; Michaelis-Menton uptake kinetics
       rnut2 = P2/(rkn(n)+P2)         !ammonium
       tnit = P1/(rkn(n)+P1)          !nitrate
       tmp = 1.0 - rnut2
!    Enforce preferential utilization of ammonium
       rnut1 = min(tmp,tnit)
       rmmn = rnut1 + rnut2
       framm = rnut2/rmmn
       rmml = tirrq/(tirrq+0.5*rikd(n))
!       rmmlice = tirrqice/(tirrqice+0.5*rikd(n))
       rmmf = P4/(rkf(n)+P4)          !iron
       rlim = min(rmml,rmmn,rmmf)
!       rlimice = min(rmmlice,rmmn,rmmf)
       grate = rmuplsr(n)*rlim*fnoice
!                + rmuplsr(n)*rlimice*(1.0-fnoice)
       rmu4(n) = grate*framm
       rmu3(n) = grate*(1.0-framm)
       gro(n) = grate*max(P(n+nnut),0.0)
       P_tend(n+nnut) = P_tend(n+nnut) + gro(n)
!    Net primary production
       pp2(n) = gro(n)*phygross*H*cchlratio*sday
#endif
#if NCHL_DEFINED > 2
!   Cyanobacteria
       n = 3
!    Nutrient-regulated growth; Michaelis-Menton uptake kinetics
       rnut2 = P2/(rkn(n)+P2)         !ammonium
       tnit = P1/(rkn(n)+P1)          !nitrate
       tmp = 1.0 - rnut2
!    Enforce preferential utilization of ammonium
       rnut1 = min(tmp,tnit)
       rmmn = rnut1 + rnut2
       framm = rnut2/rmmn
       rmml = tirrq/(tirrq+0.5*rikd(n))
       rmmf = P4/(rkf(n)+P4)          !iron
       rlim = min(rmml,rmmn,rmmf)
       rlimnfix = min(rmml,rmmf)         !limitation for N2 fixation
       rlimrkn = min(rmml,rkn(n),rmmf)   !limitation at kn
       grate = rmuplsr(n)*rlim*fnoice
       rmu4(n) = grate*framm
       rmu3(n) = grate*(1.0-framm)
       rfix = 0.25*exp(-(75.0*max(P(n+nnut),0.0)))
       rfix = max(rfix,0.0)
!       rfix = min(rfix,0.2)
       gratenfix1 = rmuplsr(n)*rlimnfix*rfix  !N fix
       graterkn = rmuplsr(n)*rlimrkn
       gratenlim = graterkn - grate
       gratenfix = min(gratenlim,gratenfix1)  !N fix cannot exceed kn
       gratenfix = max(gratenfix,0.0)*fnoice
       gron = grate*max(P(n+nnut),0.0)
       gronfix = gratenfix*max(P(n+nnut),0.0)
       gro(n) = gron + gronfix
       P_tend(n+nnut) = P_tend(n+nnut) + gro(n)
!    Net primary production
       pp2(n) = gro(n)*phygross*H*cchlratio*sday
#endif
#if NCHL_DEFINED > 3
!   Coccolithophores
       n = 4
!    Nutrient-regulated growth; Michaelis-Menton uptake kinetics
       rnut2 = P2/(rkn(n)+P2)         !ammonium
       tnit = P1/(rkn(n)+P1)          !nitrate
       tmp = 1.0 - rnut2
!    Enforce preferential utilization of ammonium
       rnut1 = min(tmp,tnit)
       rmmn = rnut1 + rnut2
       framm = rnut2/rmmn
       rmml = tirrq/(tirrq+0.5*rikd(n))
       rmmf = P4/(rkf(n)+P4)          !iron
       rlim = min(rmml,rmmn,rmmf)
       grate = rmuplsr(n)*rlim*fnoice
       rmu4(n) = grate*framm
       rmu3(n) = grate*(1.0-framm)
       gro(n) = grate*max(P(n+nnut),0.0)
       P_tend(n+nnut) = P_tend(n+nnut) + gro(n)
       gcmax = max(gcmax,grate)
!    Net primary production
       pp2(n) = gro(n)*phygross*H*cchlratio*sday
#endif
#if NCHL_DEFINED > 4
!   Dinoflagellates
       n = 5
!    Nutrient-regulated growth; Michaelis-Menton uptake kinetics
       rnut2 = P2/(rkn(n)+P2)         !ammonium
       tnit = P1/(rkn(n)+P1)          !nitrate
       tmp = 1.0 - rnut2
!    Enforce preferential utilization of ammonium
       rnut1 = min(tmp,tnit)
       rmmn = rnut1 + rnut2
       framm = rnut2/rmmn
       rmml = tirrq/(tirrq+0.5*rikd(n))
       rmmf = P4/(rkf(n)+P4)          !iron
       rlim = min(rmml,rmmn,rmmf)
       grate = rmuplsr(n)*rlim*fnoice
       rmu4(n) = grate*framm
       rmu3(n) = grate*(1.0-framm)
       gro(n) = grate*max(P(n+nnut),0.0)
       P_tend(n+nnut) = P_tend(n+nnut) + gro(n)
!    Net primary production
       pp2(n) = gro(n)*phygross*H*cchlratio*sday
#endif
#if NCHL_DEFINED > 5
!   Phaeocystis
       n = 6
!    Nutrient-regulated growth; Michaelis-Menton uptake kinetics
       rnut2 = P2/(rkn(n)+P2)         !ammonium
       tnit = P1/(rkn(n)+P1)          !nitrate
       tmp = 1.0 - rnut2
!    Enforce preferential utilization of ammonium
       rnut1 = min(tmp,tnit)
       rmmn = rnut1 + rnut2
       framm = rnut2/rmmn
       rmml = tirrq/(tirrq+0.5*rikd(n))
       rmmlice = tirrqice/(tirrqice+0.5*rikd(n))
       rmmf = P4/(rkf(n)+P4)          !iron
       rlim = min(rmml,rmmn,rmmf)
       rlimice = min(rmmlice,rmmn,rmmf)
       grate = rmuplsr(n)*rlim*fnoice                                &
                + rmuplsr(n)*rlimice*(1.0-fnoice)
       rmu4(n) = grate*framm
       rmu3(n) = grate*(1.0-framm)
       gro(n) = grate*max(P(n+nnut),0.0)
       P_tend(n+nnut) = P_tend(n+nnut) + gro(n)
!    Net primary production
       pp2(n) = gro(n)*phygross*H*cchlratio*sday
#endif
!   Nutrient uptake
       upn = 0.0
       upa = 0.0
       upf = 0.0
       do n = 1,nchl
        Pn = max(P(nnut+n),0.0)
        upn = upn + bn*(rmu3(n)*Pn)
        upa = upa + bn*(rmu4(n)*Pn)
        upf = upf + bf*gro(n)
       enddo
       P_tend(1) = P_tend(1) - upn
       P_tend(2) = P_tend(2) - upa
       ups = bs*gro(1)
       P_tend(3) = P_tend(3) - ups
       P_tend(4) = P_tend(4) - upf
      endif

      call carbon(k,cnratio,cchlratio,bn,bf,remin,P(:),tfac,          &
       fnoice,tzoo,tirrq,cdomabsq,gro(:),picdis,                      &
       atmco2,wspd,slp,T,S,H,P_tend(:),pco2,cflx)

!  Sinking rate temperature (viscosity) dependence (also convert
!   to /s)
      viscfac = 0.451 + 0.0178*T
      do n = nnut+2,npe
       ws(n) = wss(n)*viscfac*fnoice
      enddo
!      do n = nds+1,nde
!       ws(n) = wss(n)*viscfac*fnoice
!      enddo
      n = nnut+1                    !diatoms
      Pn = max(P(n),0.0)
      wsday = axs(n,2)*Pn
      wsday = min(wsday,7.0)      !block exp to < 1096
      wsday = axs(n,1)*exp(wsday)
      wsday = min(wsday,20.0)
      wss(n) = wsday/sday
      ws(n) = wss(n)*viscfac*fnoice
      n = nnut+4                    !cocco's
      ws(n) = wssc*viscfac*fnoice
      do n = nds,nde                !detrital sinking rates
       Pn = max(P(n),0.0)
       wsday = axs(n,2)*Pn
       wsday = min(wsday,7.0)      !block exp to < 1096
       wsday = axs(n,1)*exp(wsday)
       wsday = min(wsday,50.0)
       wss(n) = wsday/sday
       ws(n) = wss(n)*viscfac*fnoice
      enddo
      n = ncs+3                     !PIC
      Pn = max(P(n),0.0)
      wsday = axs(n,2)*Pn
      wsday = min(wsday,7.0)      !block exp to < 1096
      wsday = axs(n,1)*exp(wsday)
      wsday = min(wsday,50.0)
      wss(n) = wsday/sday
      ws(n) = wss(n)*viscfac*fnoice

!  Save method for hard boundary condition (no flux)
!      srate = 0.0 - wss(n)*P(i,k-1,m,n)

!  Save array for "past" layer thickness
      Hpst = H

      return
      end
