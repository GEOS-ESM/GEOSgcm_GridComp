      subroutine ppco2(slp,atmco2,T,S,P,ff,pco2)
 
!  Computes pCO2 in the surface layer and delta pCO2 with the 
!  atmosphere using OCMIP protocols.
 
#include "definebio.h"
      real, parameter :: tabar=2310.0 !mean total alkalinity uE/kg; OCMIP
      real, parameter :: stdslp=1013.25 !standard sea level pressure mb
      real, parameter :: Sbar=34.836    !global mean annual salinity 
                                        !(area-wt)
      real P(ntyp)
      data pHmin,pHmax /1.0, 16.0/
 
      dtco2 = 0.0
      pco2 = 0.0
!      atmco2 = 371.3  !uatm or ppmv (equivalent); global mean
!                      2000-2003 from OCMIP
!  Constants
      phlo = pHmin
      phhi = pHmax
      ph = 8.0
 
!  Use OCMIP subroutines
      Tsfc = T
      sal = S
      DIC = max(P(ncs+1),0.0)       !uM
      PO4 = max(P(1),0.0)*0.1       !uM phosphate, converted from NO3/PO4
!                                    ratio from Conkright et al 1994, 
!                                    global using 1st 3 standard depths 
!                                    (0, 10, and 20m)
      Si = max(P(3),0.0)            !uM Si
 
!  Convert to units for co2calc
       dic_in = DIC*1.0E-3       !uM to mol/m3
!       ta = tabar*sal/Sbar       !adjust alk for salinity
       ta = P(ncs+2)             !uM
       ta = ta/1.0245            !uM to uE/kg
       ta_in = ta*1024.5*1.0E-6  !uE/kg to E/m3
       pt_in = PO4*1.0E-3        !uM to mol/m3
       sit_in = Si*1.0E-3        !uM to mol/m3
       xco2_in = atmco2
       atmpres = slp/stdslp
 
       call co2calc_SWS(Tsfc,sal,dic_in,ta_in,pt_in,sit_in,          &
!     *            phlo,phhi,ph,xco2_in,atmpres,co2star,dco2star,
!     *            pco2surf,dpco2)                                    
                  phlo,phhi,ph,ff,pco2surf)
       pco2 = pco2surf
!       dtco2(i) = dco2star
!       pHsfc(i) = pH
 
      return
      end
