!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

      module rrtmg_sw_setcoef

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      !use parkind, only : im => kind , rb => kind 
      use parrrsw, only : mxmol
      use rrsw_ref, only : pref, preflog, tref
      use rrsw_vsn, only : hvrset, hnamset

      implicit none

      contains

!----------------------------------------------------------------------------
      subroutine setcoef_sw(ncol, nlayers, pavel, tavel, pz, tz, tbound, coldry, wkl, &
                            laytrop, layswtch, laylow, jp, jt, jt1, &
                            co2mult, colch4, colco2, colh2o, colmol, coln2o, &
                            colo2, colo3, fac00, fac01, fac10, fac11, &
                            selffac, selffrac, indself, forfac, forfrac, indfor)
!----------------------------------------------------------------------------
!
! Purpose:  For a given atmosphere, calculate the indices and
! fractions related to the pressure and temperature interpolations.

! Modifications:
! Original: J. Delamere, AER, Inc. (version 2.5, 02/04/01)
! Revised: Rewritten and adapted to ECMWF F90, JJMorcrette 030224
! Revised: For uniform rrtmg formatting, MJIacono, Jul 2006

! ------ Declarations -------

! ----- Input -----
      integer, intent(in) :: ncol

      integer , intent(in) :: nlayers         ! total number of layers
      
      real , intent(in) :: pavel(:,:)            ! layer pressures (mb) 
                                                      !    Dimensions: (nlayers)
      real , intent(in) :: tavel(:,:)            ! layer temperatures (K)
                                                      !    Dimensions: (nlayers)
      real , intent(in) :: pz(:,0:)              ! level (interface) pressures (hPa, mb)
                                                      !    Dimensions: (0:nlayers)
      real , intent(in) :: tz(:,0:)              ! level (interface) temperatures (K)
                                                      !    Dimensions: (0:nlayers)
      real , intent(in) :: tbound(:)             ! surface temperature (K)
      real , intent(in) :: coldry(:,:)           ! dry air column density (mol/cm2)
                                                      !    Dimensions: (nlayers)
      real , intent(in) :: wkl(:,:,:)            ! molecular amounts (mol/cm-2)
                                                      !    Dimensions: (mxmol,nlayers)

! ----- Output -----
      integer , intent(out) :: laytrop(:)        ! tropopause layer index
      integer , intent(out) :: layswtch(:)       ! 
      integer , intent(out) :: laylow(:)         ! 

      integer , intent(out) :: jp(:,:)           ! 
                                                      !    Dimensions: (nlayers)
      integer , intent(out) :: jt(:,:)           !
                                                      !    Dimensions: (nlayers)
      integer , intent(out) :: jt1(:,:)          !
                                                      !    Dimensions: (nlayers)

      real , intent(out) :: colh2o(:,:)          ! column amount (h2o)
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: colco2(:,:)          ! column amount (co2)
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: colo3(:,:)           ! column amount (o3)
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: coln2o(:,:)          ! column amount (n2o)
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: colch4(:,:)          ! column amount (ch4)
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: colo2(:,:)           ! column amount (o2)
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: colmol(:,:)          ! 
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: co2mult(:,:)         !
                                                      !    Dimensions: (nlayers)

      integer , intent(out) :: indself(:,:) 
                                                      !    Dimensions: (nlayers)
      integer , intent(out) :: indfor(:,:) 
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: selffac(:,:) 
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: selffrac(:,:) 
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: forfac(:,:) 
                                                      !    Dimensions: (nlayers)
      real , intent(out) :: forfrac(:,:) 
                                                      !    Dimensions: (nlayers)

      real , intent(out) :: fac00(:,:) , fac01(:,:) , fac10(:,:) , fac11(:,:)  

! ----- Local -----

      integer  :: indbound
      integer  :: indlev0
      integer  :: lay
      integer  :: jp1
      integer  :: iplon

      real  :: stpfac
      real  :: tbndfrac
      real  :: t0frac
      real  :: plog
      real  :: fp
      real  :: ft
      real  :: ft1
      real  :: water
      real  :: scalefac
      real  :: factor
      real  :: co2reg
      real  :: compfp


! Initializations
 stpfac = 296. /1013. 

   

!$acc kernels present(pavel, layswtch, laytrop, laylow)
 layswtch = 0
 laytrop = 0
 laylow = 0
 do iplon = 1, ncol
    do lay = 1, nlayers
         plog = log(pavel(iplon,lay) )
         if (plog .ge. 4.56) laytrop(iplon) = laytrop(iplon) + 1
         if (plog .ge. 6.62) laylow(iplon) = laylow(iplon) + 1
    end do
 end do
!$acc end kernels



!$acc kernels loop present(pavel, tavel, pz, tz, tbound) &
!$acc present(coldry, wkl, jp, jt, jt1, colh2o, colco2) &
!$acc present(colo3, coln2o, colch4, colo2, colmol, co2mult, indself) &
!$acc present(indfor, selffac, selffrac, forfac, forfrac, fac00, fac01, fac10, fac11)
do iplon = 1, ncol


      indbound = tbound(iplon) - 159. 
      tbndfrac = tbound(iplon) - int(tbound(iplon))
      
      
      
      indlev0  = tz(iplon,0)  - 159. 
      t0frac   = tz(iplon,0)  - int(tz(iplon,0) )


! Begin layer loop


      do lay = 1, nlayers
! Find the two reference pressures on either side of the
! layer pressure.  Store them in JP and JP1.  Store in FP the
! fraction of the difference (in ln(pressure)) between these
! two values that the layer pressure lies.

         plog = log(pavel(iplon,lay) )
         jp(iplon,lay)  = int(36.  - 5*(plog+0.04 ))
         if (jp(iplon,lay)  .lt. 1) then
            jp(iplon,lay)  = 1
         elseif (jp(iplon,lay)  .gt. 58) then
            jp(iplon,lay)  = 58
         endif
         jp1 = jp(iplon,lay)  + 1
         fp = 5.  * (preflog(jp(iplon,lay) ) - plog)

! Determine, for each reference pressure (JP and JP1), which
! reference temperature (these are different for each  
! reference pressure) is nearest the layer temperature but does
! not exceed it.  Store these indices in JT and JT1, resp.
! Store in FT (resp. FT1) the fraction of the way between JT
! (JT1) and the next highest reference temperature that the 
! layer temperature falls.

         jt(iplon,lay)  = int(3.  + (tavel(iplon,lay) -tref(jp(iplon,lay) ))/15. )
         if (jt(iplon,lay)  .lt. 1) then
            jt(iplon,lay)  = 1
         elseif (jt(iplon,lay)  .gt. 4) then
            jt(iplon,lay)  = 4
         endif
         ft = ((tavel(iplon,lay) -tref(jp(iplon,lay) ))/15. ) - float(jt(iplon,lay) -3)
         jt1(iplon,lay)  = int(3.  + (tavel(iplon,lay) -tref(jp1))/15. )
         if (jt1(iplon,lay)  .lt. 1) then
            jt1(iplon,lay)  = 1
         elseif (jt1(iplon,lay)  .gt. 4) then
            jt1(iplon,lay)  = 4
         endif
         ft1 = ((tavel(iplon,lay) -tref(jp1))/15. ) - float(jt1(iplon,lay) -3)

         water = wkl(iplon,1,lay) /coldry(iplon,lay) 
         scalefac = pavel(iplon,lay)  * stpfac / tavel(iplon,lay) 

! If the pressure is less than ~100mb, perform a different
! set of species interpolations.

         if (plog .le. 4.56 ) then

                  forfac(iplon,lay)  = scalefac / (1.+water)
         factor = (tavel(iplon,lay) -188.0 )/36.0 
         indfor(iplon,lay)  = 3
         forfrac(iplon,lay)  = factor - 1.0 

! Calculate needed column amounts.

         colh2o(iplon,lay)  = 1.e-20  * wkl(iplon,1,lay) 
         colco2(iplon,lay)  = 1.e-20  * wkl(iplon,2,lay) 
         colo3(iplon,lay)   = 1.e-20  * wkl(iplon,3,lay) 
         coln2o(iplon,lay)  = 1.e-20  * wkl(iplon,4,lay) 
         colch4(iplon,lay)  = 1.e-20  * wkl(iplon,6,lay) 
         colo2(iplon,lay)   = 1.e-20  * wkl(iplon,7,lay) 
         colmol(iplon,lay)  = 1.e-20  * coldry(iplon,lay)  + colh2o(iplon,lay) 
         if (colco2(iplon,lay)  .eq. 0. ) colco2(iplon,lay)  = 1.e-32  * coldry(iplon,lay) 
         if (coln2o(iplon,lay)  .eq. 0. ) coln2o(iplon,lay)  = 1.e-32  * coldry(iplon,lay) 
         if (colch4(iplon,lay)  .eq. 0. ) colch4(iplon,lay)  = 1.e-32  * coldry(iplon,lay) 
         if (colo2(iplon,lay)   .eq. 0. ) colo2(iplon,lay)   = 1.e-32  * coldry(iplon,lay) 
         co2reg = 3.55e-24  * coldry(iplon,lay) 
         co2mult(iplon,lay) = (colco2(iplon,lay)  - co2reg) * &
               272.63 *exp(-1919.4 /tavel(iplon,lay) )/(8.7604e-4 *tavel(iplon,lay) )

         selffac(iplon,lay)  = 0. 
         selffrac(iplon,lay) = 0. 
         indself(iplon,lay)  = 0


         else

         

! Set up factors needed to separately include the water vapor
! foreign-continuum in the calculation of absorption coefficient.

         forfac(iplon,lay)  = scalefac / (1.+water)
         factor = (332.0 -tavel(iplon,lay) )/36.0 
         indfor(iplon,lay)  = min(2, max(1, int(factor)))
         forfrac(iplon,lay)  = factor - float(indfor(iplon,lay) )

! Set up factors needed to separately include the water vapor
! self-continuum in the calculation of absorption coefficient.

         selffac(iplon,lay)  = water * forfac(iplon,lay) 
         factor = (tavel(iplon,lay) -188.0 )/7.2 
         indself(iplon,lay)  = min(9, max(1, int(factor)-7))
         selffrac(iplon,lay)  = factor - float(indself(iplon,lay)  + 7)

! Calculate needed column amounts.

         colh2o(iplon,lay)  = 1.e-20  * wkl(iplon,1,lay) 
         colco2(iplon,lay)  = 1.e-20  * wkl(iplon,2,lay) 
         colo3(iplon,lay)  = 1.e-20  * wkl(iplon,3,lay) 
!           colo3(lay) = 0. 
!           colo3(lay) = colo3(lay)/1.16 
         coln2o(iplon,lay)  = 1.e-20  * wkl(iplon,4,lay) 
         colch4(iplon,lay)  = 1.e-20  * wkl(iplon,6,lay) 
         colo2(iplon,lay)  = 1.e-20  * wkl(iplon,7,lay) 
         colmol(iplon,lay)  = 1.e-20  * coldry(iplon,lay)  + colh2o(iplon,lay) 
!           colco2(lay) = 0. 
!           colo3(lay) = 0. 
!           coln2o(lay) = 0. 
!           colch4(lay) = 0. 
!           colo2(lay) = 0. 
!           colmol(lay) = 0. 
         if (colco2(iplon,lay)  .eq. 0. ) colco2(iplon,lay)  = 1.e-32  * coldry(iplon,lay) 
         if (coln2o(iplon,lay)  .eq. 0. ) coln2o(iplon,lay)  = 1.e-32  * coldry(iplon,lay) 
         if (colch4(iplon,lay)  .eq. 0. ) colch4(iplon,lay)  = 1.e-32  * coldry(iplon,lay) 
         if (colo2(iplon,lay)  .eq. 0. ) colo2(iplon,lay)  = 1.e-32  * coldry(iplon,lay) 
! Using E = 1334.2 cm-1.
         co2reg = 3.55e-24  * coldry(iplon,lay) 
         co2mult(iplon,lay) = (colco2(iplon,lay)  - co2reg) * &
               272.63 *exp(-1919.4 /tavel(iplon,lay) )/(8.7604e-4 *tavel(iplon,lay) )
      
        end if
! We have now isolated the layer ln pressure and temperature,
! between two reference pressures and two reference temperatures 
! (for each reference pressure).  We multiply the pressure 
! fraction FP with the appropriate temperature fractions to get 
! the factors that will be needed for the interpolation that yields
! the optical depths (performed in routines TAUGBn for band n).

         compfp = 1.  - fp
         fac10(iplon,lay)  = compfp * ft
         fac00(iplon,lay)  = compfp * (1.  - ft)
         fac11(iplon,lay)  = fp * ft1
         fac01(iplon,lay)  = fp * (1.  - ft1)

! End layer loop
      end do

! End column loop
      end do
!$acc end kernels

end subroutine setcoef_sw

!***************************************************************************
      subroutine swatmref
!***************************************************************************

      save
 
! These pressures are chosen such that the ln of the first pressure
! has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
! each subsequent ln(pressure) differs from the previous one by 0.2.

      pref(:) = (/ &
          1.05363e+03 ,8.62642e+02 ,7.06272e+02 ,5.78246e+02 ,4.73428e+02 , &
          3.87610e+02 ,3.17348e+02 ,2.59823e+02 ,2.12725e+02 ,1.74164e+02 , &
          1.42594e+02 ,1.16746e+02 ,9.55835e+01 ,7.82571e+01 ,6.40715e+01 , &
          5.24573e+01 ,4.29484e+01 ,3.51632e+01 ,2.87892e+01 ,2.35706e+01 , &
          1.92980e+01 ,1.57998e+01 ,1.29358e+01 ,1.05910e+01 ,8.67114e+00 , &
          7.09933e+00 ,5.81244e+00 ,4.75882e+00 ,3.89619e+00 ,3.18993e+00 , &
          2.61170e+00 ,2.13828e+00 ,1.75067e+00 ,1.43333e+00 ,1.17351e+00 , &
          9.60789e-01 ,7.86628e-01 ,6.44036e-01 ,5.27292e-01 ,4.31710e-01 , &
          3.53455e-01 ,2.89384e-01 ,2.36928e-01 ,1.93980e-01 ,1.58817e-01 , &
          1.30029e-01 ,1.06458e-01 ,8.71608e-02 ,7.13612e-02 ,5.84256e-02 , &
          4.78349e-02 ,3.91639e-02 ,3.20647e-02 ,2.62523e-02 ,2.14936e-02 , &
          1.75975e-02 ,1.44076e-02 ,1.17959e-02 ,9.65769e-03  /)

      preflog(:) = (/ &
           6.9600e+00 , 6.7600e+00 , 6.5600e+00 , 6.3600e+00 , 6.1600e+00 , &
           5.9600e+00 , 5.7600e+00 , 5.5600e+00 , 5.3600e+00 , 5.1600e+00 , &
           4.9600e+00 , 4.7600e+00 , 4.5600e+00 , 4.3600e+00 , 4.1600e+00 , &
           3.9600e+00 , 3.7600e+00 , 3.5600e+00 , 3.3600e+00 , 3.1600e+00 , &
           2.9600e+00 , 2.7600e+00 , 2.5600e+00 , 2.3600e+00 , 2.1600e+00 , &
           1.9600e+00 , 1.7600e+00 , 1.5600e+00 , 1.3600e+00 , 1.1600e+00 , &
           9.6000e-01 , 7.6000e-01 , 5.6000e-01 , 3.6000e-01 , 1.6000e-01 , &
          -4.0000e-02 ,-2.4000e-01 ,-4.4000e-01 ,-6.4000e-01 ,-8.4000e-01 , &
          -1.0400e+00 ,-1.2400e+00 ,-1.4400e+00 ,-1.6400e+00 ,-1.8400e+00 , &
          -2.0400e+00 ,-2.2400e+00 ,-2.4400e+00 ,-2.6400e+00 ,-2.8400e+00 , &
          -3.0400e+00 ,-3.2400e+00 ,-3.4400e+00 ,-3.6400e+00 ,-3.8400e+00 , &
          -4.0400e+00 ,-4.2400e+00 ,-4.4400e+00 ,-4.6400e+00  /)

! These are the temperatures associated with the respective 
! pressures for the MLS standard atmosphere. 

      tref(:) = (/ &
           2.9420e+02 , 2.8799e+02 , 2.7894e+02 , 2.6925e+02 , 2.5983e+02 , &
           2.5017e+02 , 2.4077e+02 , 2.3179e+02 , 2.2306e+02 , 2.1578e+02 , &
           2.1570e+02 , 2.1570e+02 , 2.1570e+02 , 2.1706e+02 , 2.1858e+02 , &
           2.2018e+02 , 2.2174e+02 , 2.2328e+02 , 2.2479e+02 , 2.2655e+02 , &
           2.2834e+02 , 2.3113e+02 , 2.3401e+02 , 2.3703e+02 , 2.4022e+02 , &
           2.4371e+02 , 2.4726e+02 , 2.5085e+02 , 2.5457e+02 , 2.5832e+02 , &
           2.6216e+02 , 2.6606e+02 , 2.6999e+02 , 2.7340e+02 , 2.7536e+02 , &
           2.7568e+02 , 2.7372e+02 , 2.7163e+02 , 2.6955e+02 , 2.6593e+02 , &
           2.6211e+02 , 2.5828e+02 , 2.5360e+02 , 2.4854e+02 , 2.4348e+02 , & 
           2.3809e+02 , 2.3206e+02 , 2.2603e+02 , 2.2000e+02 , 2.1435e+02 , &
           2.0887e+02 , 2.0340e+02 , 1.9792e+02 , 1.9290e+02 , 1.8809e+02 , &
           1.8329e+02 , 1.7849e+02 , 1.7394e+02 , 1.7212e+02  /)

      end subroutine swatmref

      end module rrtmg_sw_setcoef


