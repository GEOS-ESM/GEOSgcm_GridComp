!_ ---------------------------------------------------------------------
!_ RCS lines preceded by "c_ "
!_ ---------------------------------------------------------------------
  
!_ $Source$ 
!_ $Revision$
!_ $Date$   ;  $State$
!_ $Author$ ;  $Locker$
  
!_ ---------------------------------------------------------------------
!_ $Log$
!_ Revision 1.1.2.2  2010/12/09 19:25:37  atrayano
!_ AT: added/committed files for Watson
!_
!_ Revision 1.4  2004/06/16 13:18:29  orr
!_ Initial Revision
  
!_ ---------------------------------------------------------------------
   
        subroutine ta_iter_SWS(x,fn,df)
        real k12,k12p,k123p
        real k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
        common /const/k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,hSWS
        common /species/bt,st,ft,sit,pt,dic,ta
  
! Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)
! - by A. Mouchet, 2004:
! Fixed Problems w/ version of ta_iter_1.f used in OCMIP-2 (vers. 1.2)
!  1) fixed errors in signs, parenthesis and coefficient c in derivative
!  2) changed from Total to Seawater Scale 
!     * c defined for seawater H scale; 
!     * fn and df adapted to KF on free H scale
!     * comments have been adapted
  

  
! This routine expresses TA as a function of DIC, hSWS and constants.
! It also calculates the derivative of this function with respect to 
! hSWS. It is used in the iterative solution for hSWS. In the call
! "x" is the input value for hSWS, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhSWS
  
      x2=x*x
      x3=x2*x
      k12 = k1*k2
      k12p = k1p*k2p
      k123p = k12p*k3p
      c = 1.0 + st/ks + ft/kf
      a = x3 + k1p*x2 + k12p*x + k123p
      a2=a*a
      da = 3.0*x2 + 2.0*k1p*x + k12p
      b = x2 + k1*x + k12
      b2=b*b
      db = 2.0*x + k1
  
!	fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
!===========================================================================
  
      fn = k1*x*dic/b +                                             &
             2.0*dic*k12/b +                                        &
             bt/(1.0 + x/kb) +                                      &
             kw/x +                                                 &
             pt*k12p*x/a +                                          &
             2.0*pt*k123p/a +                                       &
             sit/(1.0 + x/ksi) -                                    &
             x/c -                                                  &
             st/(1.0 + ks/(x/c)) -                                  &
             ft/(1.0 + kf/(x/c)) -                                  &
             pt*x3/a -                                              &
             ta
  
!	df = dfn/dx

!       if (a2 == 0.) then
!          a2 = tiny(a2)
!       end if

       df = ((k1*dic*b) - k1*x*dic*db)/b2 -                        &
             2.0*dic*k12*db/b2 -                                   &
             bt/kb/(1.0+x/kb)**2. -                                &
             kw/x2 +                                               &
             (pt*k12p*(a - x*da))/a2 -                             &
             2.0*pt*k123p*da/a2 -                                  &
             sit/ksi/(1.0+x/ksi)**2. -                             &
             1.0/c -                                               &
             st *(1.0 + ks/(x/c))**(-2.0) *(ks*c/x2) -             &
             ft*(1.0 + kf/(x/c))**(-2.0) *(kf*c/x2) -              &
             pt*x2*(3.0*a-x*da)/a2
  
      return
      end

