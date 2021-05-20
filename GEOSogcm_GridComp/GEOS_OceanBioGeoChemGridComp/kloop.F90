    subroutine kloop(km,DT,solFe,remin,rkn,rks,rkf,cnratio,cfratio,  &
        cchlratio,wss,wssc,drate,Rm,rlamz,greff,dratez1,               &
        dratez2,regen,axs,rmumax,                                      &
        P,H,atmfe,fnoice,                                              &
        rikd,gcmax,                                                    &
        tirrq,cdomabsq,atmco2,wspd,slp,T,S,                            &
        pco2,cflx,ppz)
!  Uncomment this to include rivers
!      subroutine kloop(km,rlat,DT,solFe,remin,rkn,rks,rkf,cnratio,cfratio,  &
!        pco2,cflx,ppz,dout,area)

!  k-loop for ocean biology

#include "definebio.h"
      real, parameter :: sday=86400.0    !seconds per day
      real rmumax(nchl)     !max phyto growth rate at 20oC, /d
      real remin(ndet)      !detrital remineralization rate /s
      real axs(ntyp,2)      !detrital sinking a and b coefficients
      real rkn(nchl),rks(nchl) !half-sat constants for N and S (uM)
      real rkf(nchl)        !half-sat constant for Fe (nM)
      real wss(ntyp)        !particle sinking rate m/s
      real P(km,ntyp)
      real gcmax(km),wssc(km)
      real H(km),tfac(km)
      real cchlratio(km)
      real atmfe,fnoice,tzoo
      real rikd(km,nchl)
      real rmuplsr(km,nchl)
      real tirrq(km),cdomabsq(km)
      real T(km),S(km)
      real P_tend(km,ntyp)
      real ws(km,ntyp)
      real pp2(nchl),ppz(nchl)
      real gronfix(km)
      real zpic(km),picdis(km)
      data csd /3500.0/  !Ca compensation depth from OCMIP
      data disrat /0.17/ !dissolution rate d-1 from Buitenhuis etal 2001
      real dout,area
      real P_tendr(km,ntyp)

!  Bio Loop
      gronfix = 0.0
      P_tend = 0.0
      P_tendr = 0.0
      ppz = 0.0
!      ylat = rlat*57.29578
      k = 1
      do while (k .le. km )
       if (H(k) >= 1.0E10) exit
       kend = k
       k = k+1
      enddo
!  Dissolution of PIC
      picdis = 0.0
      k = 1
      zpic = 0.0
      zc = 75.0       !from Buitenhuis et al 2001
      zpic(k) = H(k)
      do k = 2,kend
       zpic(k) = zpic(k-1) + H(k)
       if (zpic(k) .gt. csd)then
        picdis(k) = disrat/sday*exp((zpic(k)-zc)/csd)
       endif
      enddo
!  Loop
      do k = 1,kend
       call ptend(k,DT,solFe,remin,rkn,rks,rkf,cnratio,cfratio,        &
        cchlratio(k),wss,wssc(k),drate,Rm,rlamz,greff,dratez1,dratez2, &
        regen,axs,                                                     &
        P(k,:),H(k),atmfe,fnoice,rikd(k,:),rmumax,gcmax(k),            &
        tirrq(k),cdomabsq(k),atmco2,wspd,slp,T(k),S(k),                &
        P_tend(k,:),gronfix(k),ws(k,:),Hpst,picdis(k),pco2,cflx,pp2)
!   Restore N from N-fixation into lower layers (reverse from uptake)
       kto = (kend-k)+1
       bn = cchlratio(k)/cnratio
       P_tend(kto,1:1) = P_tend(kto,1:1) - bn*gronfix(k)
       ppz = ppz + pp2
       fnoice = 1.0    !no ice processes below surface
      enddo

      call sink(km,kend,H,ws,P,cnratio,P_tend)
      P = P + P_tend*DT

!! Add nutrient in river discharge
!       if (dout>0.0 .and. dout<0.1E16) then
!         call rivers(km,ylat,DT,P,H,dout,area,P_tendr)
!         P(1,1:ntyp) = P(1,1:ntyp) + P_tendr(1,1:ntyp)*DT
!!         write(6,'(A,E12.5,E12.5,E12.5)') &
!!         " Updated nitrate in kloop is ", P(1,1),P_tend(1,1),P_tendr(1,1)
!       endif

      return
      end
