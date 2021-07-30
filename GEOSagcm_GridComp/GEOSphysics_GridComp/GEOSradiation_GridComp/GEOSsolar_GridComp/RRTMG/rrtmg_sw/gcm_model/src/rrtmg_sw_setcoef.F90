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

   use parrrsw, only : mxmol
   use rrsw_ref, only : pref, preflog, tref

   implicit none

contains

   !----------------------------------------------------------------------------
   subroutine setcoef_sw( &
      ncol, nlayers, pavel, tavel, coldry, wkl, &
      laytrop, laylow, jp, jt, jt1, &
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
   !
   !----------------------------------------------------------------------------

      ! ----- Input -----

      integer, intent(in) :: ncol            ! number of gridcols
      integer, intent(in) :: nlayers         ! number of layers
      
      real, intent(in) :: pavel(:,:)            ! layer pressures (mb) 
                                                      !    Dimensions: (nlayers)
      real, intent(in) :: tavel(:,:)            ! layer temperatures (K)
                                                      !    Dimensions: (nlayers)
      real, intent(in) :: coldry(:,:)           ! dry air column density (mol/cm2)
                                                      !    Dimensions: (nlayers)
      real, intent(in) :: wkl(:,:,:)            ! molecular amounts (mol/cm-2)
                                                      !    Dimensions: (mxmol,nlayers)
      ! ----- Output -----

      integer, intent(out) :: laytrop(:)        ! tropopause layer index
      integer, intent(out) :: laylow(:)         ! 

      integer, intent(out) :: jp(:,:)           ! 
                                                      !    Dimensions: (nlayers)
      integer, intent(out) :: jt(:,:)           !
                                                      !    Dimensions: (nlayers)
      integer, intent(out) :: jt1(:,:)          !
                                                      !    Dimensions: (nlayers)

      real, intent(out) :: colh2o(:,:)          ! column amount (h2o)
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: colco2(:,:)          ! column amount (co2)
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: colo3(:,:)           ! column amount (o3)
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: coln2o(:,:)          ! column amount (n2o)
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: colch4(:,:)          ! column amount (ch4)
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: colo2(:,:)           ! column amount (o2)
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: colmol(:,:)          ! 
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: co2mult(:,:)         !
                                                      !    Dimensions: (nlayers)

      integer, intent(out) :: indself(:,:) 
                                                      !    Dimensions: (nlayers)
      integer, intent(out) :: indfor(:,:) 
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: selffac(:,:) 
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: selffrac(:,:) 
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: forfac(:,:) 
                                                      !    Dimensions: (nlayers)
      real, intent(out) :: forfrac(:,:) 
                                                      !    Dimensions: (nlayers)

      real, intent(out) :: fac00(:,:), fac01(:,:), fac10(:,:), fac11(:,:)  

      ! ----- Local -----

      integer :: icol, lay, jp1
      real :: plog, fp, ft, ft1, water, scalefac, factor, co2reg, compfp

      ! Initializations
      real, parameter :: stpfac = 296. / 1013. 

      laytrop = 0
      laylow = 0
      do icol = 1,ncol
         do lay = 1,nlayers
            plog = log(pavel(icol,lay))
            if (plog >= 4.56) laytrop(icol) = laytrop(icol) + 1
            if (plog >= 6.62) laylow(icol) = laylow(icol) + 1
         end do
      end do

      do icol = 1, ncol
         do lay = 1, nlayers

            ! Find the two reference pressures on either side of the
            ! layer pressure.  Store them in JP and JP1. Store in FP the
            ! fraction of the difference (in ln(pressure)) between these
            ! two values that the layer pressure lies.

            plog = log(pavel(icol,lay))
            jp(icol,lay) = int(36. - 5*(plog+0.04 ))
            if (jp(icol,lay) < 1) then
               jp(icol,lay) = 1
            elseif (jp(icol,lay) > 58) then
               jp(icol,lay) = 58
            endif
            jp1 = jp(icol,lay) + 1
            fp = 5. * (preflog(jp(icol,lay) ) - plog)

            ! Determine, for each reference pressure (JP and JP1), which
            ! reference temperature (these are different for each  
            ! reference pressure) is nearest the layer temperature but does
            ! not exceed it. Store these indices in JT and JT1, resp.
            ! Store in FT (resp. FT1) the fraction of the way between JT
            ! (JT1) and the next highest reference temperature that the 
            ! layer temperature falls.

            jt(icol,lay) = int(3. + (tavel(icol,lay)-tref(jp(icol,lay)))/15.)
            if (jt(icol,lay) < 1) then
               jt(icol,lay) = 1
            elseif (jt(icol,lay) > 4) then
               jt(icol,lay) = 4
            endif
            ft = ((tavel(icol,lay)-tref(jp(icol,lay)))/15.) - float(jt(icol,lay)-3)
            jt1(icol,lay) = int(3. + (tavel(icol,lay)-tref(jp1))/15.)
            if (jt1(icol,lay) < 1) then
               jt1(icol,lay) = 1
            elseif (jt1(icol,lay) > 4) then
               jt1(icol,lay) = 4
            endif
            ft1 = ((tavel(icol,lay)-tref(jp1))/15.) - float(jt1(icol,lay)-3)

            water = wkl(icol,1,lay) / coldry(icol,lay) 
            scalefac = pavel(icol,lay)  * stpfac / tavel(icol,lay) 

            ! If the pressure is less than ~100mb, perform a different
            ! set of species interpolations.

            if (plog <= 4.56 ) then

               forfac(icol,lay) = scalefac / (1.+water)
               factor = (tavel(icol,lay)-188.) / 36. 
               indfor(icol,lay) = 3
               forfrac(icol,lay) = factor - 1. 

               ! Calculate needed column amounts.

               colh2o(icol,lay) = 1.e-20 * wkl(icol,1,lay) 
               colco2(icol,lay) = 1.e-20 * wkl(icol,2,lay) 
               colo3 (icol,lay) = 1.e-20 * wkl(icol,3,lay) 
               coln2o(icol,lay) = 1.e-20 * wkl(icol,4,lay) 
               colch4(icol,lay) = 1.e-20 * wkl(icol,6,lay) 
               colo2 (icol,lay) = 1.e-20 * wkl(icol,7,lay) 
               colmol(icol,lay) = 1.e-20 * coldry(icol,lay) + colh2o(icol,lay) 
               if (colco2(icol,lay) == 0.) colco2(icol,lay) = 1.e-32 * coldry(icol,lay) 
               if (coln2o(icol,lay) == 0.) coln2o(icol,lay) = 1.e-32 * coldry(icol,lay) 
               if (colch4(icol,lay) == 0.) colch4(icol,lay) = 1.e-32 * coldry(icol,lay) 
               if (colo2 (icol,lay) == 0.) colo2 (icol,lay) = 1.e-32 * coldry(icol,lay) 
               co2reg = 3.55e-24 * coldry(icol,lay) 
               co2mult(icol,lay) = (colco2(icol,lay) - co2reg) * &
                  272.63 * exp(-1919.4/tavel(icol,lay)) / (8.7604e-4 * tavel(icol,lay))

               selffac(icol,lay)  = 0. 
               selffrac(icol,lay) = 0. 
               indself(icol,lay)  = 0

            else

               ! Set up factors needed to separately include the water vapor
               ! foreign-continuum in the calculation of absorption coefficient.

               forfac(icol,lay) = scalefac / (1.+water)
               factor = (332.-tavel(icol,lay))/36. 
               indfor(icol,lay) = min(2,max(1,int(factor)))
               forfrac(icol,lay) = factor - float(indfor(icol,lay))

               ! Set up factors needed to separately include the water vapor
               ! self-continuum in the calculation of absorption coefficient.

               selffac(icol,lay) = water * forfac(icol,lay) 
               factor = (tavel(icol,lay)-188.)/7.2 
               indself(icol,lay) = min(9,max(1,int(factor)-7))
               selffrac(icol,lay) = factor - float(indself(icol,lay) + 7)

               ! Calculate needed column amounts.

               colh2o(icol,lay) = 1.e-20 * wkl(icol,1,lay) 
               colco2(icol,lay) = 1.e-20 * wkl(icol,2,lay) 
               colo3 (icol,lay) = 1.e-20 * wkl(icol,3,lay) 
               coln2o(icol,lay) = 1.e-20 * wkl(icol,4,lay) 
               colch4(icol,lay) = 1.e-20 * wkl(icol,6,lay) 
               colo2 (icol,lay) = 1.e-20 * wkl(icol,7,lay) 
               colmol(icol,lay) = 1.e-20 * coldry(icol,lay) + colh2o(icol,lay) 
               if (colco2(icol,lay) == 0.) colco2(icol,lay) = 1.e-32 * coldry(icol,lay) 
               if (coln2o(icol,lay) == 0.) coln2o(icol,lay) = 1.e-32 * coldry(icol,lay) 
               if (colch4(icol,lay) == 0.) colch4(icol,lay) = 1.e-32 * coldry(icol,lay) 
               if (colo2 (icol,lay) == 0.) colo2 (icol,lay) = 1.e-32 * coldry(icol,lay) 
               ! Using E = 1334.2 cm-1.
               co2reg = 3.55e-24 * coldry(icol,lay) 
               co2mult(icol,lay) = (colco2(icol,lay) - co2reg) * &
                  272.63 * exp(-1919.4/tavel(icol,lay))/(8.7604e-4 * tavel(icol,lay))
      
            end if

            ! We have now isolated the layer ln pressure and temperature,
            ! between two reference pressures and two reference temperatures 
            ! (for each reference pressure). We multiply the pressure 
            ! fraction FP with the appropriate temperature fractions to get 
            ! the factors that will be needed for the interpolation that yields
            ! the optical depths (performed in routines TAUGBn for band n).

            compfp = 1. - fp
            fac10(icol,lay) = compfp * ft
            fac00(icol,lay) = compfp * (1.-ft)
            fac11(icol,lay) = fp * ft1
            fac01(icol,lay) = fp * (1.-ft1)

         end do  ! layer loop
      end do  ! column loop

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


