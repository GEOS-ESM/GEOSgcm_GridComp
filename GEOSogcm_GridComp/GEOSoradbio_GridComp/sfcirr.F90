      subroutine sfcirr(lam,Fobar,thray,oza,awv,ao,aco2,           &
       asl,bsl,csl,dsl,esl,fsl,ica,                                &
       daycor,cosunz,pres,ws,ozone,wvapor,relhum,                  &
       ta,wa,asym,am,Vi,cov,cldtau,clwp,re,Ed,Es)

!  Calls clrtrans.F to get cloud-free transmittance and slingo.F to
!  get cloudy transmittance, then computes total irradiance in
!  W/m2/(variable)nm weighted by the cloudiness.

!  Tdclr is spectral clear sky direct transmittance
!  Tsclr is spectral clear sky diffuse transmittance
!  Tdcld is spectral cloudy direct transmittance
!  Tscld is spectral cloudy diffuse transmittance

#include "comlte.h"
      save ifst,p0,ozfac1,ozfac2
      integer lam(nlt)
      real Fobar(nlt),thray(nlt),oza(nlt),awv(nlt),ao(nlt),aco2(nlt)
      real asl(ncld),bsl(ncld),csl(ncld),dsl(ncld),esl(ncld),fsl(ncld)
      integer ica(nlt)
      real ta(nlt),wa(nlt),asym(nlt),Tgas(nlt)
      real Tdclr(nlt),Tsclr(nlt)
      real Tdcld(ncld),Tscld(ncld)
      real Ed(nlt),Es(nlt)
      real Edclr(nlt),Esclr(nlt)
      real Edcld(nlt),Escld(nlt)
      data ifst /0/

!  Initial constants
      if (ifst .eq. 0)then
       ozfac1 = 44.0/6370.0
       ozfac2 = 1.0 + 22.0/6370.0
       p0 = 1013.25
       ifst = 1
      endif

      Ed = 0.0
      Es = 0.0
      if (pres .lt. 0.0 .or. ws .lt. 0.0 .or. relhum .lt. 0.0     &
      .or. ozone .lt. 0.0 .or. wvapor .lt. 0.0)then
       go to 5
      endif

!  Compute atmospheric path lengths (air mass); not pressure-corrected
      sunz = acos(cosunz)*57.29578
      rtmp = (93.885-sunz)**(-1.253)
      rmu0 = cosunz+0.15*rtmp
      rm = 1.0/rmu0
      otmp = (cosunz*cosunz+ozfac1)**0.5
      rmo = ozfac2/otmp

!  Compute pressure-corrected atmospheric path length (air mass)
      rmp = pres/p0*rm

!  Loop to compute total irradiances at each grid point
!   Compute direct and diffuse irradiance for a cloudy and non-cloudy
!   atmosphere
!   Account for gaseous absorption
      do nl = 1,nlt
!    Ozone
       to = oza(nl)*ozone*0.001   !optical thickness
       oarg = -to*rmo
!   Oxygen/gases
       ag = ao(nl) + aco2(nl)
       gtmp = (1.0 + 118.3*ag*rmp)**0.45
       gtmp2 = -1.41*ag*rmp
       garg = gtmp2/gtmp
!   Water Vapor
       wtmp = (1.0+20.07*awv(nl)*wvapor*rm)**0.45
       wtmp2 = -0.2385*awv(nl)*wvapor*rm
       warg = wtmp2/wtmp
       Tgas(nl) = exp(oarg+garg+warg)
      enddo

!  Compute clear sky transmittances
      call clrtrans(lam,thray,cosunz,rm,rmp,rmo,ws,relhum,am,Vi,    &
                    ta,wa,asym,Tdclr,Tsclr)
      do nl = 1,nlt
       Foinc = Fobar(nl)*daycor*cosunz
!    Direct irradiance 
       Edir = Foinc*Tgas(nl)*Tdclr(nl)
!    Diffuse irradiance
       Edif = Foinc*Tgas(nl)*Tsclr(nl)
!    Spectral components
       Edclr(nl) = Edir
       Esclr(nl) = Edif
      enddo    !end clear nl loop
!  Compute cloudy transmittances
!      call slingo(rmu0,cldtc(ic,jc),rlwp(ic,jc),cdre(ic,jc),
!     *            Tdcld,Tscld)
      call slingo(asl,bsl,csl,dsl,esl,fsl,                            &
       rmu0,cldtau,clwp,re,Tdcld,Tscld)
      do nl = 1,nlt
       Foinc = Fobar(nl)*daycor*cosunz
!    Direct irradiance 
       Edir = Foinc*Tgas(nl)*Tdcld(ica(nl))
!    Diffuse irradiance
       Edif = Foinc*Tgas(nl)*Tscld(ica(nl))
!    Spectral components
       Edcld(nl) = Edir
       Escld(nl) = Edif
      enddo   !end cloudy nl loop

!  Sum clear and cloudy by percent
      ccov1 = cov*0.01  !convert from percent to fraction
      do nl = 1,nlt
       Ed(nl) = (1.0-ccov1)*Edclr(nl) + ccov1*Edcld(nl)
       Es(nl) = (1.0-ccov1)*Esclr(nl) + ccov1*Escld(nl)
      enddo

!   Total short-wave (W/m2)
!       sirrsw = sirrsw + Edir + Edif

5     continue
      return
      end
