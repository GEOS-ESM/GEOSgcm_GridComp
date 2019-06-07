      subroutine slingo(asl,bsl,csl,dsl,esl,fsl,                      &
                        rmu0,cldtau,clwp,cre,Tcd,Tcs)

!  Slingo's (1989) Delta-Eddington approximation for the two-
!  stream equations applied to clouds.  
!  Inputs:
!       rmu0    Kasten's approx for cosine of solar zenith angle
!       cldtau  cloud optical thickness (at 0.6 um)
!       clwp    liquid water path in cloud (g/m2)
!       cre     cloud droplet effective radius (um)
!  Outputs
!       ica     index for translating cloud arrays to light arrays
!       Tcd     spectral transmittance for downwelling direct irradiance
!       Tcs     spectral transmittance for downwelling diffuse irradiance
!  Internal
!       Tcu     spectral transmittance for upwelling irradiance

!  This version uses global mean re for the last 2 bands 2.91-3.42 and
!  3.42-4.00 um to account for unknown errors in the GEOS-5/ESMF system


#include "comlte.h"
      real asl(ncld),bsl(ncld),csl(ncld),dsl(ncld)
      real esl(ncld),fsl(ncld)
      real Tcd(ncld),Tcs(ncld)
      real Tcu(ncld)

      Tcd = 0.0
      Tcs = 0.0
      Tcu = 0.0
      
!  Compute re as funtion of cldtau and LWP according to eq. 1 in 
!  Slingo.
!   tau is derived at this wavelength (0.6 um) in the ISCCP data set
!      re = clwp*bsl(9)/(cldtau - clwp*asl(9))
!      re = min(re,15.0)  !block high re -- produces excessive direct
!  Changes to the ISCCP-D2 data set make this relationship untenable
!  (excessive re's are derived).  Instead choose a fixed re of 10 um
!  for ocean (Kiehl et al., 1998 -- J. Clim.)
!       re = 10.0
!  Paper by Han et al., 1994 (J.Clim.) show mean ocean cloud radius
!  = 11.8 um
!       re = 11.8
!  Mean of Kiehl and Han
      re = (10.0+11.8)/2.0
      remean = re

!  Compute spectral cloud characteristics
!   If MODIS re is available use it; otherwise use parameterized re above
      if (cre .ge. 0.0)then   !use MODIS re
       re = cre
      endif
      do nc = 1,22
       tauc = clwp*(asl(nc)+bsl(nc)/re)
       oneomega = csl(nc) + dsl(nc)*re
       omega = 1.0 - oneomega
       g = esl(nc) + fsl(nc)*re
       call slingomath(tauc,omega,g,rmu0,Tcd(nc),Tcs(nc))
      enddo

!  Slingo bands 23 and 24 fail due to re.  Workaround is to use mean
!  ocean re
      do nc = 23,24
       tauc = clwp*(asl(nc)+bsl(nc)/remean)
       oneomega = csl(nc) + dsl(nc)*remean
       omega = 1.0 - oneomega
       g = esl(nc) + fsl(nc)*remean
       call slingomath(tauc,omega,g,rmu0,Tcd(nc),Tcs(nc))
      enddo

      return
      end

!******************************************************************
      subroutine slingomath(tauc,omega,g,rmu0,Tcds,Tcss)

!  Mathematical calculations for the Slingo model

      real*8 alpha1,alpha2,eps,sqarg
      real*8 E,rM,Esq,rMsq,EM,val1,val2,val3

      U1 = 7.0/4.0

      b0 = 3.0/7.0*(1.0-g)
      bmu0 = 0.5 - 0.75*rmu0*g/(1.0+g)
      f = g*g
      U2 = U1*(1.0-((1.0-omega)/(7.0*omega*b0)))
      U2 = max(U2,0.0)
      alpha1 = U1*(1.0D0-omega*(1.0D0-b0))
      alpha2 = U2*omega*b0
      alpha3 = (1.0-f)*omega*bmu0
      alpha4 = (1.0-f)*omega*(1.0-bmu0)
      sqarg = alpha1*alpha1 - alpha2*alpha2
      sqarg = max(sqarg,1.0D-17)
      eps = sqrt(sqarg)
      rM = alpha2/(alpha1+eps)
      E = exp(-eps*tauc)
      val1 = 1.0D0 - omega*f
      val2 = eps*eps*rmu0*rmu0
      rnum = val1*alpha3 - rmu0*(alpha1*alpha3+alpha2*alpha4)
      rden = val1*val1 - val2
      gama1 = rnum/rden
      rnum = -val1*alpha4 - rmu0*(alpha1*alpha4+alpha2*alpha3)
      gama2 = rnum/rden
      Tdb = exp(-val1*tauc/rmu0)
      Esq = E*E
      rMsq = rM*rM
      EM = Esq*rMsq
      val3 = 1.0D0 - EM
      Rdif = rM*(1.0-Esq)/val3
      Tdif = E*(1.0-rMsq)/val3
!      Rdir = -gama2*Rdif - gama1*Tdb*Tdif + gama1
      Tdir = -gama2*Tdif - gama1*Tdb*Rdif + gama2*Tdb
!      Tdir = max(Tdir,0.0)
      Tcds = Tdb
      Tcss = max(Tdir,0.0)
!       Tcu = Tdif

      return
      end
