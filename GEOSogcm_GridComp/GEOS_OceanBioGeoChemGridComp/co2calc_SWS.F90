!_ ---------------------------------------------------------------------
!_ RCS lines preceded by "c_ "
!_ --------------------------------------------------------------------
!_
!_ $Source$ 
!_ $Revision$
!_ $Date$   ;  $State$
!_ $Author$ ;  $Locker$
!_
!_ ---------------------------------------------------------------------
!_ $Log$
!_ Revision 1.1.2.2  2010/12/09 19:25:37  atrayano
!_ AT: added/committed files for Watson
!_
!_ Revision 1.5  2004/07/02 15:29:01  orr
!_ Added missing parenthesis at end of "kf" equation (bug identified by Anne Mouchet)
!_
!_ Revision 1.4  2004/06/16 13:18:29  orr
!_ Initial Revision
!_
!_ ---------------------------------------------------------------------
!_ 
      subroutine co2calc_SWS(t,s,dic_in,ta_in,pt_in,sit_in            &
!     &                  ,phlo,phhi,ph,xco2_in,atmpres
!     &                  ,co2star,dco2star,pCO2surf,dpco2)
                        ,phlo,phhi,ph,ffco2,pCO2surf)
 
!-------------------------------------------------------------------------
 
! Modified from co2calc.f (RCS version 1.8, OCMIP-2) 
! - by A. Mouchet, 2004:
 
! NOW; "All" constants are given on seawater H scale (hSWS) 
! - HOWEVER, dissociation constants for S and F are on the 'free' H scale 
!            (this is necessary to work with hSWS)
! - this routine corrects for inconsistencies in the previous version.
 
! - Other changes:
!   * use ta_iter_SWS instead of ta_iter_1;
!   * hSWS replaces htotal;
!   * moved upward the calculation of bt, st and ft 
!     (needed in evaluation of kb);
!   * added various comments and references.
 
 
! SUBROUTINE CO2CALC_SWS
 
! PURPOSE
!	Calculate delta co2* from total alkalinity and total CO2 at
! temperature (t), salinity (s) and "atmpres" atmosphere total pressure. 
 
! USAGE
!       call co2calc_SWS(t,s,dic_in,ta_in,pt_in,sit_in
!    &                  ,phlo,phhi,ph,xco2_in,atmpres
!    &                  ,co2star,dco2star,pCO2surf,dpco2)
 
! INPUT
!	dic_in = total inorganic carbon (mol/m^3) 
!                where 1 T = 1 metric ton = 1000 kg
!	ta_in  = total alkalinity (eq/m^3) 
!	pt_in  = inorganic phosphate (mol/m^3) 
!	sit_in = inorganic silicate (mol/m^3) 
!	t      = temperature (degrees C)
!	s      = salinity (PSU)
!	phlo   = lower limit of pH range
!	phhi   = upper limit of pH range
!	xco2_in=atmospheric mole fraction CO2 in dry air (ppmv) 
!	atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
 
!       Note: arguments dic_in, ta_in, pt_in, sit_in, and xco2_in are 
!             used to initialize variables dic, ta, pt, sit, and xco2.
!             * Variables dic, ta, pt, and sit are in the common block 
!               "species".
!             * Variable xco2 is a local variable.
!             * Variables with "_in" suffix have different units 
!               than those without.

! OUTPUT
!	co2star  = CO2*water (mol/m^3)
!	dco2star = delta CO2 (mol/m^3)
!       pco2surf = oceanic pCO2 (ppmv)
!       dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
 
! IMPORTANT: Some words about units - (JCO, 4/4/1999)
!     - Models carry tracers in mol/m^3 (on a per volume basis)
!     - Conversely, this routine, which was written by observationalists 
!       (C. Sabine and R. Key), passes input arguments in umol/kg  
!       (i.e., on a per mass basis)
!     - I have changed things slightly so that input arguments are in mol/m^3,
!     - Thus, all input concentrations (dic_in, ta_in, pt_in, and st_in) 
!       should be given in mol/m^3; output arguments "co2star" and "dco2star"  
!       are likewise in mol/m^3.

! FILES and PROGRAMS NEEDED
!	drtsafe
!	ta_iter_SWS
 
!--------------------------------------------------------------------------
 
        real invtk,is,is2
        real k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
        common /const/k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,hSWS
        common /species/bt,st,ft,sit,pt,dic,ta
        external ta_iter_SWS
 

!       ---------------------------------------------------------------------
!       Change units from the input of mol/m^3 -> mol/kg:
!       (1 mol/m^3)  x (1 m^3/1024.5 kg)
!       where the ocean's mean surface density is 1024.5 kg/m^3
!       Note: mol/kg are actually what the body of this routine uses 
!       for calculations.  
!       ---------------------------------------------------------------------
        permil = 1.0 / 1024.5
!       To convert input in mol/m^3 -> mol/kg 
        pt=pt_in*permil
        sit=sit_in*permil
        ta=ta_in*permil
        dic=dic_in*permil

!       ---------------------------------------------------------------------
!       Change units from uatm to atm. That is, atm is what the body of 
!       this routine uses for calculations.
!       ---------------------------------------------------------------------
        permeg=1.e-6
!!       To convert input in uatm -> atm    !WG
!        xco2=xco2_in*permeg    !WG
!       ---------------------------------------------------------------------
 
!*************************************************************************
! Calculate all constants needed to convert between various measured
! carbon species. References for each equation are noted in the code. 
! Once calculated, the constants are
! stored and passed in the common block "const". The original version of this
! code was based on the code by Dickson in Version 2 of "Handbook of Methods
! for the Analysis of the Various Parameters of the Carbon Dioxide System
! in Seawater", DOE, 1994 (SOP No. 3, p25-26). 
 
! Derive simple terms used more than once
 
      tk = 273.15 + t
      tk100 = tk/100.0
      tk1002=tk100*tk100
      invtk=1.0/tk
      dlogtk=log(tk)
      is=19.924*s/(1000.-1.005*s)
      is2=is*is
      sqrtis=sqrt(is)
      s2=s*s
      sqrts=sqrt(s)
      s15=s**1.5
      scl=s/1.80655
 
!------------------------------------------------------------------------
! Calculate concentrations for borate, sulfate, and fluoride
 
! Uppstrom (1974)
      bt = 0.000232 * scl/10.811
! Morris & Riley (1966)
      st = 0.14 * scl/96.062
! Riley (1965)
      ft = 0.000067 * scl/18.9984
 
!------------------------------------------------------------------------
! f = k0(1-pH2O)*correction term for non-ideality
 
! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
 
      ff = exp(-162.8301 + 218.2968/tk100  +                         &
       90.9241*log(tk100) - 1.47696*tk1002 +                         &
       s * (.025695 - .025225*tk100 +                                &
       0.0049867*tk1002))
 
! K0 from Weiss 1974
 
      k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) +      &
       s * (.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))

 
!------------------------------------------------------------------------
! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3]     on hSWS
 
! Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale 
! (Original reference: Dickson and Millero, DSR, 1987)
 
      k1=10**(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk -            &
       0.0118 * s + 0.000116*s2))
 
      k2=10**(-1*(1394.7*invtk + 4.777 -                             &
       0.0184*s + 0.000118*s2))
 
!------------------------------------------------------------------------
! k1p = [H][H2PO4]/[H3PO4] on hSWS
 
! Millero p.670 (1995)
 
      k1p = exp(-4576.752*invtk + 115.540 - 18.453 * dlogtk +       &
       (-106.736*invtk + 0.69171) * sqrts +                         &
       (-0.65643*invtk - 0.01844) * s)
 
!------------------------------------------------------------------------
! k2p = [H][HPO4]/[H2PO4] on hSWS
 
! Millero p.670 (1995)
 
      k2p = exp(-8814.715*invtk + 172.1033 - 27.927 * dlogtk +      &
       (-160.340*invtk + 1.3566) * sqrts +                          &
       (0.37335*invtk - 0.05778) * s)
 
!------------------------------------------------------------------------
! k3p = [H][PO4]/[HPO4] on hSWS
 
! Millero p.670 (1995)
 
      k3p = exp(-3070.75*invtk - 18.126 +                           &
       (17.27039*invtk + 2.81197) *                                 &
       sqrts + (-44.99486*invtk - 0.09984) * s)
 
!------------------------------------------------------------------------
! ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
 
! Millero p.671 (1995) using data from Yao and Millero (1995)
! change to (mol/ kg soln)
 
      ksi = exp(-8904.2*invtk + 117.400 - 19.334 * dlogtk +          &
       (-458.79*invtk + 3.5913) * sqrtis +                           &
       (188.74*invtk - 1.5998) * is +                                &
       (-12.1652*invtk + 0.07871) * is2 +                            &
       log(1.0-0.001005*s))
 
!------------------------------------------------------------------------
! kw = [H][OH] on hSWS
 
! Millero p.670 (1995) using composite data
 
      kw = exp(-13847.26*invtk + 148.9802 - 23.6521 * dlogtk +      &
       (118.67*invtk - 5.977 + 1.0495 * dlogtk) *                   &
       sqrts - 0.01615 * s)
 
!------------------------------------------------------------------------
! ks = [H][SO4]/[HSO4] on free H scale
 
! Dickson (1990, J. chem. Thermodynamics 22, 113)
! change to (mol/ kg soln)
 
      ks=exp(-4276.1*invtk + 141.328 - 23.093*dlogtk +               &
       (-13856*invtk + 324.57 - 47.986*dlogtk) * sqrtis +            &
       (35474*invtk - 771.54 + 114.723*dlogtk) * is -                &
       2698*invtk*is**1.5 + 1776*invtk*is2 +                         &
       log(1.0 - 0.001005*s))
 
!------------------------------------------------------------------------
! kf = [H][F]/[HF] on free H scale
 
! Dickson and Riley (1979)
! change to (mol/ kg soln)
 
      kf=exp(1590.2*invtk - 12.641 + 1.525*sqrtis +                 &
       log(1.0 - 0.001005*s)) 
 
!------------------------------------------------------------------------
! kb = [H][BO2]/[HBO2] on hSWS
 
! Dickson p.673 (1990)
! change from htotal to hSWS
 
      kb=exp( (-8966.90 - 2890.53*sqrts - 77.942*s +               &
       1.728*s15 - 0.0996*s2)*invtk +                              &
       (148.0248 + 137.1942*sqrts + 1.62142*s) +                   &
       (-24.4344 - 25.085*sqrts - 0.2474*s) *                      &
       dlogtk + 0.053105*sqrts*tk +                                &
       log((1+(st/ks)+(ft/kf))/(1+(st/ks))) )
 
!*************************************************************************
 
! Calculate [H+] SWS when DIC and TA are known at T, S and 1 atm.
! The solution converges to err of xacc. The solution must be within
! the range x1 to x2.
 
! If DIC and TA are known then either a root finding or iterative method
! must be used to calculate hSWS. In this case we use the Newton-Raphson
! "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
! error trapping removed).
 
! As currently set, this procedure iterates about 12 times. The x1 and x2
! values set below will accomodate ANY oceanographic values. If an initial
! guess of the pH is known, then the number of iterations can be reduced to
! about 5 by narrowing the gap between x1 and x2. It is recommended that
! the first few time steps be run with x1 and x2 set as below. After that,
! set x1 and x2 to the previous value of the pH +/- ~0.5. The current
! setting of xacc will result in co2star accurate to 3 significant figures
! (xx.y). Making xacc bigger will result in faster convergence also, but this
! is not recommended (xacc of 10**-9 drops precision to 2 significant figures).
 
! Parentheses added around negative exponents (Keith Lindsay)
 
      x1 = 10.0**(-phhi)
      x2 = 10.0**(-phlo)
!	xacc = 1.e-10
      xacc = 1.e-14
      hSWS = drtsafe(ta_iter_SWS,x1,x2,xacc)
 
! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
 
      hSWS2=hSWS*hSWS
      co2star=dic*hSWS2/(hSWS2 + k1*hSWS + k1*k2)
!      co2starair=xco2*ff*atmpres    !WG
!      dco2star=co2starair-co2star    !WG
      ph=-log10(hSWS)

 
!       ---------------------------------------------------------------
!!      Add two output arguments for storing pCO2surf
!!      Should we be using K0 or ff for the solubility here?
!       ---------------------------------------------------------------
        pCO2surf = co2star / ff
        ffco2 = ff
!        dpCO2    = pCO2surf - xco2*atmpres    !WG
 
!  Convert units of output arguments
!      Note: co2star and dco2star are calculated in mol/kg within this routine 
!      Thus Convert now from mol/kg -> mol/m^3
!       co2star  = co2star / permil    !WG
!       dco2star = dco2star / permil    !WG

!      Note: pCO2surf and dpCO2 are calculated in atm above. 
!      Thus convert now to uatm
       pCO2surf = pCO2surf / permeg
!       dpCO2    = dpCO2 / permeg    !WG
 
      return
      end
