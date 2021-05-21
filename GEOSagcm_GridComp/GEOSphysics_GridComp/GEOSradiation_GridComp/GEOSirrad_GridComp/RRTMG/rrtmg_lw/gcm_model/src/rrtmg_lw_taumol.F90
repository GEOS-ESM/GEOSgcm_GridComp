module rrtmg_lw_taumol

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

   use parrrtm,  only : nbndlw, ngptlw
   use rrlw_wvn, only : nspa, nspb, ngb
   use rrlw_con, only : oneminus
   use rrtmg_lw_setcoef

   implicit none

! *******************************************************************************
! *                                                                             *
! *                  Optical depths developed for the                           *
! *                                                                             *
! *                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
! *                                                                             *
! *                                                                             *
! *            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
! *                        131 HARTWELL AVENUE                                  *
! *                        LEXINGTON, MA 02421                                  *
! *                                                                             *
! *                                                                             *
! *                           ELI J. MLAWER                                     * 
! *                         JENNIFER DELAMERE                                   * 
! *                         STEVEN J. TAUBMAN                                   *
! *                         SHEPARD A. CLOUGH                                   *
! *                                                                             *
! *                                                                             *
! *                                                                             *
! *                                                                             *
! *                       email:  mlawer@aer.com                                *
! *                       email:  jdelamer@aer.com                              *
! *                                                                             *
! *        The authors wish to acknowledge the contributions of the             *
! *        following people:  Karen Cady-Pereira, Patrick D. Brown,             *  
! *        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom.    *
! *                                                                             *
! *******************************************************************************
! *                                                                             *
! *  Revision for g-point reduction: Michael J. Iacono, AER, Inc.               *
! *                                                                             *
! *******************************************************************************
! *     TAUMOL                                                                  *
! *                                                                             *
! *     This file contains the subroutines TAUGBn (where n goes from            *
! *     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
! *     per g-value and layer for band n.                                       *
! *                                                                             *
! *  Output:  optical depths (unitless)                                         *
! *           fractions needed to compute Planck functions at every layer       *
! *               and g-value                                                   *
! *                                                                             *
! *     COMMON /TAUGCOM/  taug(MXLAY,MG)                                        *
! *     COMMON /PLANKG/   pfracs(MXLAY,MG)                                      *
! *                                                                             *
! *  Input                                                                      *
! *                                                                             *
! *     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
! *     COMMON /PRECISE/  oneminus                                              *
! *     COMMON /PROFILE/  nlay,PAVEL(MXLAY),TAVEL(MXLAY),                       *
! *     &                 PZ(0:MXLAY),TZ(0:MXLAY)                               *
! *     COMMON /PROFDATA/ LAYTROP,                                              *
! *    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),             *
! *    &                  COLN2O(MXLAY),colco(MXLAY),COLCH4(MXLAY),             *
! *    &                  COLO2(MXLAY)                                          *
! *     COMMON /INTFAC/   fac00(icol,MXLAY),fac01(icol,MXLAY),                  *
! *    &                  FAC10(MXLAY),fac11(icol,MXLAY)                        *
! *     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
! *     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
! *                                                                             *
! *     Description:                                                            *
! *     NG(IBAND) - number of g-values in band IBAND                            *
! *     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
! *                   atmospheres that are stored for band IBAND per            *
! *                   pressure level and temperature.  Each of these            *
! *                   atmospheres has different relative amounts of the         *
! *                   key species for the band (i.e. different binary           *
! *                   species parameters).                                      *
! *     NSPB(IBAND) - same for upper atmosphere                                 *
! *     oneminus - since problems are caused in some cases by interpolation     *
! *                parameters equal to or greater than 1, for these cases       *
! *                these parameters are set to this value, slightly < 1.        *
! *     PAVEL - layer pressures (mb)                                            *
! *     TAVEL - layer temperatures (degrees K)                                  *
! *     PZ - level pressures (mb)                                               *
! *     TZ - level temperatures (degrees K)                                     *
! *     LAYTROP - layer at which switch is made from one combination of         *
! *               key species to another                                        *
! *     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
! *               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
! *               respectively (molecules/cm**2)                                *
! *     FACij(LAY) - for layer LAY, these are factors that are needed to        *
! *                  compute the interpolation factors that multiply the        *
! *                  appropriate reference k-values.  A value of 0 (1) for      *
! *                  i,j indicates that the corresponding factor multiplies     *
! *                  reference k-value for the lower (higher) of the two        *
! *                  appropriate temperatures, and altitudes, respectively.     *
! *     JP - the index of the lower (in altitude) of the two appropriate        *
! *          reference pressure levels needed for interpolation                 *
! *     JT, JT1 - the indices of the lower of the two appropriate reference     *
! *               temperatures needed for interpolation (for pressure           *
! *               levels JP and JP+1, respectively)                             *
! *     SELFFAC - scale factor needed for water vapor self-continuum, equals    *
! *               (water vapor density)/(atmospheric density at 296K and        *
! *               1013 mb)                                                      *
! *     SELFFRAC - factor needed for temperature interpolation of reference     *
! *                water vapor self-continuum data                              *
! *     INDSELF - index of the lower of the two appropriate reference           *
! *               temperatures needed for the self-continuum interpolation      *
! *     FORFAC  - scale factor needed for water vapor foreign-continuum.        *
! *     FORFRAC - factor needed for temperature interpolation of reference      *
! *                water vapor foreign-continuum data                           *
! *     INDFOR  - index of the lower of the two appropriate reference           *
! *               temperatures needed for the foreign-continuum interpolation   *
! *                                                                             *
! *  Data input                                                                 *
! *     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG),*
! *                 FORREF(4,MG), KA_M'MGAS', KB_M'MGAS'                        *
! *        (note:  n is the band number,'MGAS' is the species name of the minor *
! *         gas)                                                                *
! *                                                                             *
! *     Description:                                                            *
! *     KA - k-values for low reference atmospheres (key-species only)          *
! *          (units: cm**2/molecule)                                            *
! *     KB - k-values for high reference atmospheres (key-species only)         *
! *          (units: cm**2/molecule)                                            *
! *     KA_M'MGAS' - k-values for low reference atmosphere minor species        *
! *          (units: cm**2/molecule)                                            *
! *     KB_M'MGAS' - k-values for high reference atmosphere minor species       *
! *          (units: cm**2/molecule)                                            *
! *     SELFREF - k-values for water vapor self-continuum for reference         *
! *               atmospheres (used below LAYTROP)                              *
! *               (units: cm**2/molecule)                                       *
! *     FORREF  - k-values for water vapor foreign-continuum for reference      *
! *               atmospheres (used below/above LAYTROP)                        *
! *               (units: cm**2/molecule)                                       *
! *                                                                             *
! *     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
! *     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
! *                                                                             *
!********************************************************************************

contains

   !--------------------------------------------------------
   subroutine taumol (ncol, nlay, pavel, taua, taug, pfracs)
   !--------------------------------------------------------

      integer, intent(in)  :: ncol                       ! number of columns
      integer, intent(in)  :: nlay                       ! number of layers
      real,    intent(in)  :: pavel  (nlay,ncol)         ! layer pressures [hPa]
      real,    intent(in)  :: taua   (nlay,nbndlw,ncol)  ! aerosol optical depths
      real,    intent(out) :: taug   (nlay,ngptlw,ncol)  ! gas optical depths
      real,    intent(out) :: pfracs (nlay,ngptlw,ncol)  ! Planck fractions

      ! Calculate gaseous optical depth and Planck fractions
      !   for each spectral band
      call taugb1  (ncol, nlay, taug, pfracs, pavel)
      call taugb2  (ncol, nlay, taug, pfracs, pavel)
      call taugb3  (ncol, nlay, taug, pfracs)
      call taugb4  (ncol, nlay, taug, pfracs)
      call taugb5  (ncol, nlay, taug, pfracs)
      call taugb6  (ncol, nlay, taug, pfracs)
      call taugb7  (ncol, nlay, taug, pfracs)
      call taugb8  (ncol, nlay, taug, pfracs)
      call taugb9  (ncol, nlay, taug, pfracs)
      call taugb10 (ncol, nlay, taug, pfracs)
      call taugb11 (ncol, nlay, taug, pfracs)
      call taugb12 (ncol, nlay, taug, pfracs)
      call taugb13 (ncol, nlay, taug, pfracs)
      call taugb14 (ncol, nlay, taug, pfracs)
      call taugb15 (ncol, nlay, taug, pfracs)
      call taugb16 (ncol, nlay, taug, pfracs)

      ! add the aerosol optical depths to the gas optical depths
      call addAerosols (ncol, nlay, taua, taug)

   end subroutine taumol


   !--------------------------------------------------
   subroutine taugb1 (ncol, nlay, taug, pfracs, pavel)
   !--------------------------------------------------

   !  Written by Eli J. Mlawer, Atmospheric & Environmental Research.
   !  Revised by Michael J. Iacono, Atmospheric & Environmental Research.
   !
   !     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
   !                          (high key - h2o; high minor - n2)
   !
   !        note: previous versions of rrtm band 1: 
   !           10-250 cm-1 (low - h2o; high - h2o)

   ! Compute the optical depth by interpolating in ln(pressure) and 
   ! temperature. Below laytrop, the water vapor self-continuum and
   ! foreign continuum is interpolated (in temperature) separately.

      use rrlw_kg01

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)
      real,    intent(in)    :: pavel(nlay,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, indm, ig
      real :: pp, corradj, scalen2, tauself, taufor, taun2
     
      ! Minor gas mapping levels:
      !     lower - n2, p = 142.5490 mb, t = 215.70 K
      !     upper - n2, p = 142.5490 mb, t = 215.70 K

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(1) + 1
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(1) + 1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)
         pp = pavel(lay,icol)
         corradj =  1.
         if (pp .lt. 250.) then
            corradj = 1. - 0.15 * (250. - pp) / 154.4 
         endif
         scalen2 = colbrd(lay,icol) * scaleminorn2(lay,icol)

         do ig = 1,ng1
            tauself = selffac(lay,icol) * (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) -  forref(indf,ig))) 
            taun2 = scalen2 * (ka_mn2(indm,ig) + & 
               minorfrac(lay,icol) * (ka_mn2(indm+1,ig) - ka_mn2(indm,ig)))
            taug(lay,ig,icol) = corradj * (colh2o(lay,icol) * &
               (fac00(lay,icol) * absa(ind0,  ig) + &
                fac10(lay,icol) * absa(ind0+1,ig) + &
                fac01(lay,icol) * absa(ind1,  ig) + &
                fac11(lay,icol) * absa(ind1+1,ig))  & 
               + tauself + taufor + taun2)
            pfracs(lay,ig,icol) = fracrefa(ig)
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(1) + 1
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(1) + 1
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)
         pp = pavel(lay,icol)
         corradj = 1. - 0.15 * (pp / 95.6)
         scalen2 = colbrd(lay,icol) * scaleminorn2(lay,icol)

         do ig = 1,ng1
            taufor = forfac(lay,icol) * (forref(indf,ig) + &
               forfrac(lay,icol) * (forref(indf+1,ig) - forref(indf,ig))) 
            taun2 = scalen2 * (kb_mn2(indm,ig) + & 
               minorfrac(lay,icol) * (kb_mn2(indm+1,ig) - kb_mn2(indm,ig)))
            taug(lay,ig,icol) = corradj * (colh2o(lay,icol) * &
               (fac00(lay,icol) * absb(ind0,  ig) + &
                fac10(lay,icol) * absb(ind0+1,ig) + &
                fac01(lay,icol) * absb(ind1,  ig) + &
                fac11(lay,icol) * absb(ind1+1,ig))  &  
               + taufor + taun2)
            pfracs(lay,ig,icol) = fracrefb(ig)
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb1


   !--------------------------------------------------
   subroutine taugb2 (ncol, nlay, taug, pfracs, pavel)
   !--------------------------------------------------
   !
   !     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
   !
   !     note: previous version of rrtm band 2: 
   !           250 - 500 cm-1 (low - h2o; high - h2o)

   ! Compute the optical depth by interpolating in ln(pressure) and 
   ! temperature. Below laytrop, the water vapor self-continuum and
   ! foreign continuum is interpolated (in temperature) separately.

      use parrrtm, only : ngs1
      use rrlw_kg02

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)
      real,    intent(in)    :: pavel(nlay,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, ig
      real :: pp, corradj, tauself, taufor

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(2) + 1
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(2) + 1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         pp = pavel(lay,icol)
         corradj = 1. - .05 * (pp - 100.) / 900. 

         do ig = 1,ng2
            tauself = selffac(lay,icol) * (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            taug(lay,ngs1+ig,icol) = corradj * (colh2o(lay,icol) * &
               (fac00(lay,icol) * absa(ind0,  ig) + &
                fac10(lay,icol) * absa(ind0+1,ig) + &
                fac01(lay,icol) * absa(ind1,  ig) + &
                fac11(lay,icol) * absa(ind1+1,ig))  &
               + tauself + taufor)
            pfracs(lay,ngs1+ig,icol) = fracrefa(ig)
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(2) + 1
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(2) + 1
         indf = indfor(lay,icol)

         do ig = 1,ng2
            taufor = forfac(lay,icol) * (forref(indf,ig) + &
               forfrac(lay,icol) * (forref(indf+1,ig) - forref(indf,ig))) 
            taug(lay,ngs1+ig,icol) = colh2o(lay,icol) * &
               (fac00(lay,icol) * absb(ind0,  ig) + &
                fac10(lay,icol) * absb(ind0+1,ig) + &
                fac01(lay,icol) * absb(ind1,  ig) + &
                fac11(lay,icol) * absb(ind1+1,ig))  &
               + taufor
            pfracs(lay,ngs1+ig,icol) = fracrefb(ig)
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb2


   !-------------------------------------------
   subroutine taugb3 (ncol, nlay, taug, pfracs)
   !-------------------------------------------
   !
   !     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
   !                           (high key - h2o,co2; high minor - n2o)

   ! Compute the optical depth by interpolating in ln(pressure) and 
   ! temperature, and appropriate species. Below laytrop, the water vapor 
   ! self-continuum and foreign continuum is interpolated (in temperature) 
   ! separately.

      use parrrtm, only : ngs2
      use rrlw_ref, only : chi_mls
      use rrlw_kg03

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmn2o, jpl
      real :: speccomb, specparm, specmult, fs
      real :: speccomb1, specparm1, specmult1, fs1
      real :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, &
              fmn2o, fmn2omf, chi_n2o, ratn2o, adjfac, adjcoln2o
      real :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real :: p, p4, fk0, fk1, fk2
      real :: fac000, fac100, fac200, fac010, fac110, fac210
      real :: fac001, fac101, fac201, fac011, fac111, fac211
      real :: tauself, taufor, n2om1, n2om2, absn2o
      real :: refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b
      real :: tau_major, tau_major1

!     pref(:) = (/ &
!         1.05363e+03 ,8.62642e+02 ,7.06272e+02 ,5.78246e+02 ,4.73428e+02 , &
!         3.87610e+02 ,3.17348e+02 ,2.59823e+02 ,2.12725e+02 ,1.74164e+02 , &
!         1.42594e+02 ,1.16746e+02 ,9.55835e+01 ,7.82571e+01 ,6.40715e+01 , &
!         5.24573e+01 ,4.29484e+01 ,3.51632e+01 ,2.87892e+01 ,2.35706e+01 , &
!         1.92980e+01 ,1.57998e+01 ,1.29358e+01 ,1.05910e+01 ,8.67114e+00 , &
!         7.09933e+00 ,5.81244e+00 ,4.75882e+00 ,3.89619e+00 ,3.18993e+00 , &
!         2.61170e+00 ,2.13828e+00 ,1.75067e+00 ,1.43333e+00 ,1.17351e+00 , &
!         9.60789e-01 ,7.86628e-01 ,6.44036e-01 ,5.27292e-01 ,4.31710e-01 , &
!         3.53455e-01 ,2.89384e-01 ,2.36928e-01 ,1.93980e-01 ,1.58817e-01 , &
!         1.30029e-01 ,1.06458e-01 ,8.71608e-02 ,7.13612e-02 ,5.84256e-02 , &
!         4.78349e-02 ,3.91639e-02 ,3.20647e-02 ,2.62523e-02 ,2.14936e-02 , &
!         1.75975e-02 ,1.44076e-02 ,1.17959e-02 ,9.65769e-03 /)

      ! Minor gas mapping levels:
      !     lower - n2o, p = 706.272 mb, t = 278.94 K
      !     upper - n2o, p = 95.58 mb, t = 215.7 K

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower/upper atmosphere.

      ! P = 212.725 mb (level 9)
      refrat_planck_a = chi_mls(1,9)/chi_mls(2,9)

      ! P = 706.270 mb (level 3)
      refrat_m_a = chi_mls(1,3)/chi_mls(2,3)

      ! P = 95.58 mb
      refrat_planck_b = chi_mls(1,13)/chi_mls(2,13)
      refrat_m_b = refrat_planck_b

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         speccomb = colh2o(lay,icol) + rat_h2oco2(lay,icol)*colco2(lay,icol)
         specparm = colh2o(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )        

         speccomb1 = colh2o(lay,icol) + rat_h2oco2_1(lay,icol)*colco2(lay,icol)
         specparm1 = colh2o(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mn2o = colh2o(lay,icol) + refrat_m_a*colco2(lay,icol)
         specparm_mn2o = colh2o(lay,icol)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminus) specparm_mn2o = oneminus
         specmult_mn2o = 8. *specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0 )
         fmn2omf = minorfrac(lay,icol)*fmn2o
!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(lay,icol)/coldry(lay,icol)
         ratn2o = 1.e20 *chi_n2o/chi_mls(4,jp(lay,icol)+1)
         if (ratn2o .gt. 1.5 ) then
            adjfac = 0.5 +(ratn2o-0.5 )**0.65 
            adjcoln2o = adjfac*chi_mls(4,jp(lay,icol)+1)*coldry(lay,icol)*1.e-20 
         else
            adjcoln2o = coln2o(lay,icol)
         endif

         speccomb_planck = colh2o(lay,icol)+refrat_planck_a*colco2(lay,icol)
         specparm_planck = colh2o(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(3) + js
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(3) + js1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else
            fac000 = (1.  - fs) * fac00(lay,icol)
            fac010 = (1.  - fs) * fac10(lay,icol)
            fac100 = fs * fac00(lay,icol)
            fac110 = fs * fac10(lay,icol)
         endif
         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else
            fac001 = (1.  - fs1) * fac01(lay,icol)
            fac011 = (1.  - fs1) * fac11(lay,icol)
            fac101 = fs1 * fac01(lay,icol)
            fac111 = fs1 * fac11(lay,icol)
         endif

         do ig = 1,ng3
            tauself = selffac(lay,icol)* (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
               (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
            n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
               (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(lay,icol) * (n2om2 - n2om1)

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac200 * absa(ind0+2, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig) + &
                   fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                  (fac200 * absa(ind0-1, ig) + &
                   fac100 * absa(ind0,   ig) + &
                   fac000 * absa(ind0+1, ig) + &
                   fac210 * absa(ind0+8, ig) + &
                   fac110 * absa(ind0+9, ig) + &
                   fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac201 * absa(ind1+2, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig) + &
                   fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                  (fac201 * absa(ind1-1, ig) + &
                   fac101 * absa(ind1,   ig) + &
                   fac001 * absa(ind1+1, ig) + &
                   fac211 * absa(ind1+8, ig) + &
                   fac111 * absa(ind1+9, ig) + &
                   fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs2+ig,icol) = tau_major + tau_major1 &
               + tauself + taufor + adjcoln2o * absn2o
            pfracs(lay,ngs2+ig,icol) = fracrefa(ig,jpl) + fpl * &
               (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
    
      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         speccomb = colh2o(lay,icol) + rat_h2oco2(lay,icol)*colco2(lay,icol)
         specparm = colh2o(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(lay,icol) + rat_h2oco2_1(lay,icol)*colco2(lay,icol)
         specparm1 = colh2o(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 4. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         fac000 = (1.  - fs) * fac00(lay,icol)
         fac010 = (1.  - fs) * fac10(lay,icol)
         fac100 = fs * fac00(lay,icol)
         fac110 = fs * fac10(lay,icol)
         fac001 = (1.  - fs1) * fac01(lay,icol)
         fac011 = (1.  - fs1) * fac11(lay,icol)
         fac101 = fs1 * fac01(lay,icol)
         fac111 = fs1 * fac11(lay,icol)

         speccomb_mn2o = colh2o(lay,icol) + refrat_m_b*colco2(lay,icol)
         specparm_mn2o = colh2o(lay,icol)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminus) specparm_mn2o = oneminus
         specmult_mn2o = 4. *specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0 )
         fmn2omf = minorfrac(lay,icol)*fmn2o
!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(lay,icol)/coldry(lay,icol)
         ratn2o = 1.e20*chi_n2o/chi_mls(4,jp(lay,icol)+1)
         if (ratn2o .gt. 1.5 ) then
            adjfac = 0.5 +(ratn2o-0.5 )**0.65 
            adjcoln2o = adjfac*chi_mls(4,jp(lay,icol)+1)*coldry(lay,icol)*1.e-20 
         else
            adjcoln2o = coln2o(lay,icol)
         endif

         speccomb_planck = colh2o(lay,icol)+refrat_planck_b*colco2(lay,icol)
         specparm_planck = colh2o(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 4. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(3) + js
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(3) + js1
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)

         do ig = 1,ng3
            taufor = forfac(lay,icol) * (forref(indf,ig) + &
               forfrac(lay,icol) * (forref(indf+1,ig) - forref(indf,ig))) 
            n2om1 = kb_mn2o(jmn2o,indm,ig) + fmn2o * &
               (kb_mn2o(jmn2o+1,indm,ig)-kb_mn2o(jmn2o,indm,ig))
            n2om2 = kb_mn2o(jmn2o,indm+1,ig) + fmn2o * &
               (kb_mn2o(jmn2o+1,indm+1,ig)-kb_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(lay,icol) * (n2om2 - n2om1)
            taug(lay,ngs2+ig,icol) = &
               speccomb * &
               (fac000 * absb(ind0,  ig) + &
                fac100 * absb(ind0+1,ig) + &
                fac010 * absb(ind0+5,ig) + &
                fac110 * absb(ind0+6,ig))  &
               + speccomb1 * &
               (fac001 * absb(ind1,  ig) + &
                fac101 * absb(ind1+1,ig) + &
                fac011 * absb(ind1+5,ig) + &
                fac111 * absb(ind1+6,ig))  &
               + taufor + adjcoln2o * absn2o
            pfracs(lay,ngs2+ig,icol) = fracrefb(ig,jpl) + fpl * &
               (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb3


   !-------------------------------------------
   subroutine taugb4 (ncol, nlay, taug, pfracs)
   !-------------------------------------------
   !
   !     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)

   ! Compute the optical depth by interpolating in ln(pressure) and 
   ! temperature, and appropriate species. Below laytrop, the water 
   ! vapor self-continuum and foreign continuum is interpolated (in temperature) 
   ! separately.

      use parrrtm, only : ngs3
      use rrlw_ref, only : chi_mls
      use rrlw_kg04

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, ig
      integer :: js, js1, jpl
      real :: speccomb, specparm, specmult, fs
      real :: speccomb1, specparm1, specmult1, fs1
      real :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real :: p, p4, fk0, fk1, fk2
      real :: fac000, fac100, fac200, fac010, fac110, fac210
      real :: fac001, fac101, fac201, fac011, fac111, fac211
      real :: tauself, taufor
      real :: refrat_planck_a, refrat_planck_b
      real :: tau_major, tau_major1

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower/upper atmosphere.

      ! P = 142.5940 mb
      refrat_planck_a = chi_mls(1,11)/chi_mls(2,11)

      ! P = 95.58350 mb
      refrat_planck_b = chi_mls(3,13)/chi_mls(2,13)

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         speccomb = colh2o(lay,icol) + rat_h2oco2(lay,icol)*colco2(lay,icol)
         specparm = colh2o(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(lay,icol) + rat_h2oco2_1(lay,icol)*colco2(lay,icol)
         specparm1 = colh2o(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_planck = colh2o(lay,icol)+refrat_planck_a*colco2(lay,icol)
         specparm_planck = colh2o(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(4) + js
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(4) + js1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else
            fac000 = (1.  - fs) * fac00(lay,icol)
            fac010 = (1.  - fs) * fac10(lay,icol)
            fac100 = fs * fac00(lay,icol)
            fac110 = fs * fac10(lay,icol)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else
            fac001 = (1.  - fs1) * fac01(lay,icol)
            fac011 = (1.  - fs1) * fac11(lay,icol)
            fac101 = fs1 * fac01(lay,icol)
            fac111 = fs1 * fac11(lay,icol)
         endif

         do ig = 1,ng4
            tauself = selffac(lay,icol)* (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac200 * absa(ind0+2, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig) + &
                   fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                  (fac200 * absa(ind0-1, ig) + &
                   fac100 * absa(ind0,   ig) + &
                   fac000 * absa(ind0+1, ig) + &
                   fac210 * absa(ind0+8, ig) + &
                   fac110 * absa(ind0+9, ig) + &
                   fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac201 * absa(ind1+2, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig) + &
                   fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                  (fac201 * absa(ind1-1, ig) + &
                   fac101 * absa(ind1,   ig) + &
                   fac001 * absa(ind1+1, ig) + &
                   fac211 * absa(ind1+8, ig) + &
                   fac111 * absa(ind1+9, ig) + &
                   fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs3+ig,icol) = &
               tau_major + tau_major1 + tauself + taufor
            pfracs(lay,ngs3+ig,icol) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
    
      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         speccomb = colo3(lay,icol) + rat_o3co2(lay,icol)*colco2(lay,icol)
         specparm = colo3(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colo3(lay,icol) + rat_o3co2_1(lay,icol)*colco2(lay,icol)
         specparm1 = colo3(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 4. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         fac000 = (1.  - fs) * fac00(lay,icol)
         fac010 = (1.  - fs) * fac10(lay,icol)
         fac100 = fs * fac00(lay,icol)
         fac110 = fs * fac10(lay,icol)
         fac001 = (1.  - fs1) * fac01(lay,icol)
         fac011 = (1.  - fs1) * fac11(lay,icol)
         fac101 = fs1 * fac01(lay,icol)
         fac111 = fs1 * fac11(lay,icol)

         speccomb_planck = colo3(lay,icol)+refrat_planck_b*colco2(lay,icol)
         specparm_planck = colo3(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 4. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(4) + js
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(4) + js1

         do ig = 1,ng4
            taug(lay,ngs3+ig,icol) = &
               speccomb * &
               (fac000 * absb(ind0,  ig) + &
                fac100 * absb(ind0+1,ig) + &
                fac010 * absb(ind0+5,ig) + &
                fac110 * absb(ind0+6,ig))  &
               + speccomb1 * &
               (fac001 * absb(ind1,  ig) + &
                fac101 * absb(ind1+1,ig) + &
                fac011 * absb(ind1+5,ig) + &
                fac111 * absb(ind1+6,ig))
            pfracs(lay,ngs3+ig,icol) = fracrefb(ig,jpl) + fpl * &
               (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
         enddo

! Empirical modification to code to improve stratospheric cooling rates
! for co2. Revised to apply weighting for g-point reduction in this band.

         taug(lay,ngs3+8,icol)=taug(lay,ngs3+8,icol)*0.92
         taug(lay,ngs3+9,icol)=taug(lay,ngs3+9,icol)*0.88
         taug(lay,ngs3+10,icol)=taug(lay,ngs3+10,icol)*1.07
         taug(lay,ngs3+11,icol)=taug(lay,ngs3+11,icol)*1.1
         taug(lay,ngs3+12,icol)=taug(lay,ngs3+12,icol)*0.99
         taug(lay,ngs3+13,icol)=taug(lay,ngs3+13,icol)*0.88
         taug(lay,ngs3+14,icol)=taug(lay,ngs3+14,icol)*0.943

      end do  ! upper layers

      end do  ! column

   end subroutine taugb4


   !-------------------------------------------
   subroutine taugb5 (ncol, nlay, taug, pfracs)
   !-------------------------------------------
   !
   !     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
   !                           (high key - o3,co2)

   ! Compute the optical depth by interpolating in ln(pressure) and 
   ! temperature, and appropriate species. Below laytrop, the 
   ! water vapor self-continuum and foreign continuum is 
   ! interpolated (in temperature) separately.

      use parrrtm, only : ngs4
      use rrlw_ref, only : chi_mls
      use rrlw_kg05

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmo3, jpl
      real :: speccomb, specparm, specmult, fs
      real :: speccomb1, specparm1, specmult1, fs1
      real :: speccomb_mo3, specparm_mo3, specmult_mo3, fmo3
      real :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real :: p, p4, fk0, fk1, fk2
      real :: fac000, fac100, fac200, fac010, fac110, fac210
      real :: fac001, fac101, fac201, fac011, fac111, fac211
      real :: tauself, taufor, o3m1, o3m2, abso3
      real :: refrat_planck_a, refrat_planck_b, refrat_m_a
      real :: tau_major, tau_major1

      ! Minor gas mapping levels:
      !     lower - o3, p = 317.34 mb, t = 240.77 K
      !     lower - ccl4

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower/upper atmosphere.

      ! P = 473.420 mb
      refrat_planck_a = chi_mls(1,5)/chi_mls(2,5)

      ! P = 317.3480
      refrat_m_a = chi_mls(1,7)/chi_mls(2,7)

      ! P = 0.2369 mb
      refrat_planck_b = chi_mls(3,43)/chi_mls(2,43)

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         speccomb = colh2o(lay,icol) + rat_h2oco2(lay,icol)*colco2(lay,icol)
         specparm = colh2o(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(lay,icol) + rat_h2oco2_1(lay,icol)*colco2(lay,icol)
         specparm1 = colh2o(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mo3 = colh2o(lay,icol) + refrat_m_a*colco2(lay,icol)
         specparm_mo3 = colh2o(lay,icol)/speccomb_mo3
         if (specparm_mo3 .ge. oneminus) specparm_mo3 = oneminus
         specmult_mo3 = 8. *specparm_mo3
         jmo3 = 1 + int(specmult_mo3)
         fmo3 = mod(specmult_mo3,1.0 )

         speccomb_planck = colh2o(lay,icol)+refrat_planck_a*colco2(lay,icol)
         specparm_planck = colh2o(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(5) + js
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(5) + js1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else
            fac000 = (1.  - fs) * fac00(lay,icol)
            fac010 = (1.  - fs) * fac10(lay,icol)
            fac100 = fs * fac00(lay,icol)
            fac110 = fs * fac10(lay,icol)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else
            fac001 = (1.  - fs1) * fac01(lay,icol)
            fac011 = (1.  - fs1) * fac11(lay,icol)
            fac101 = fs1 * fac01(lay,icol)
            fac111 = fs1 * fac11(lay,icol)
         endif

         do ig = 1,ng5
            tauself = selffac(lay,icol) * (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            o3m1 = ka_mo3(jmo3,indm,ig) + fmo3 * &
               (ka_mo3(jmo3+1,indm,ig)-ka_mo3(jmo3,indm,ig))
            o3m2 = ka_mo3(jmo3,indm+1,ig) + fmo3 * &
               (ka_mo3(jmo3+1,indm+1,ig)-ka_mo3(jmo3,indm+1,ig))
            abso3 = o3m1 + minorfrac(lay,icol)*(o3m2-o3m1)

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac200 * absa(ind0+2, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig) + &
                   fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                  (fac200 * absa(ind0-1, ig) + &
                   fac100 * absa(ind0,   ig) + &
                   fac000 * absa(ind0+1, ig) + &
                   fac210 * absa(ind0+8, ig) + &
                   fac110 * absa(ind0+9, ig) + &
                   fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac201 * absa(ind1+2, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig) + &
                   fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * & 
                  (fac201 * absa(ind1-1, ig) + &
                   fac101 * absa(ind1,   ig) + &
                   fac001 * absa(ind1+1, ig) + &
                   fac211 * absa(ind1+8, ig) + &
                   fac111 * absa(ind1+9, ig) + &
                   fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs4+ig,icol) = &
               tau_major + tau_major1 + tauself + taufor &
                 + abso3 * colo3(lay,icol) &
                 + colccl4(lay,icol) * ccl4(ig)
            pfracs(lay,ngs4+ig,icol) = fracrefa(ig,jpl) + fpl * &
               (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         speccomb = colo3(lay,icol) + rat_o3co2(lay,icol)*colco2(lay,icol)
         specparm = colo3(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colo3(lay,icol) + rat_o3co2_1(lay,icol)*colco2(lay,icol)
         specparm1 = colo3(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 4. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         fac000 = (1.  - fs) * fac00(lay,icol)
         fac010 = (1.  - fs) * fac10(lay,icol)
         fac100 = fs * fac00(lay,icol)
         fac110 = fs * fac10(lay,icol)
         fac001 = (1.  - fs1) * fac01(lay,icol)
         fac011 = (1.  - fs1) * fac11(lay,icol)
         fac101 = fs1 * fac01(lay,icol)
         fac111 = fs1 * fac11(lay,icol)

         speccomb_planck = colo3(lay,icol)+refrat_planck_b*colco2(lay,icol)
         specparm_planck = colo3(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 4. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(5) + js
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(5) + js1
         
         do ig = 1,ng5
            taug(lay,ngs4+ig,icol) = &
               speccomb * &
               (fac000 * absb(ind0,  ig) + &
                fac100 * absb(ind0+1,ig) + &
                fac010 * absb(ind0+5,ig) + &
                fac110 * absb(ind0+6,ig)) &
               + speccomb1 * &
               (fac001 * absb(ind1,  ig) + &
                fac101 * absb(ind1+1,ig) + &
                fac011 * absb(ind1+5,ig) + &
                fac111 * absb(ind1+6,ig))  &
               + colccl4(lay,icol) * ccl4(ig)
            pfracs(lay,ngs4+ig,icol) = fracrefb(ig,jpl) + fpl * &
               (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb5


   !-------------------------------------------
   subroutine taugb6 (ncol, nlay, taug, pfracs)
   !-------------------------------------------
   !
   !     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
   !                           (high key - nothing; high minor - cfc11, cfc12)

   ! Compute the optical depth by interpolating in ln(pressure) and
   ! temperature. The water vapor self-continuum and foreign continuum
   ! is interpolated (in temperature) separately.  

      use parrrtm, only : ngs5
      use rrlw_ref, only : chi_mls
      use rrlw_kg06

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, indm, ig
      real :: chi_co2, ratco2, adjfac, adjcolco2
      real :: tauself, taufor, absco2
     
      ! Minor gas mapping levels:
      !     lower - co2, p = 706.2720 mb, t = 294.2 K
      !     upper - cfc11, cfc12

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

! In atmospheres where the amount of CO2 is too great to be considered
! a minor species, adjust the column amount of CO2 by an empirical factor 
! to obtain the proper contribution.
         chi_co2 = colco2(lay,icol)/(coldry(lay,icol))
         ratco2 = 1.e20 *chi_co2/chi_mls(2,jp(lay,icol)+1)
         if (ratco2 .gt. 3.0 ) then
            adjfac = 2.0 +(ratco2-2.0 )**0.77 
            adjcolco2 = adjfac*chi_mls(2,jp(lay,icol)+1)*coldry(lay,icol)*1.e-20 
         else
            adjcolco2 = colco2(lay,icol)
         endif

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(6) + 1
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(6) + 1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)

         do ig = 1,ng6
            tauself = selffac(lay,icol) * (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig)))
            absco2 = (ka_mco2(indm,ig) + minorfrac(lay,icol) * &
               (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            taug(lay,ngs5+ig,icol) = colh2o(lay,icol) * &
               (fac00(lay,icol) * absa(ind0,  ig) + &
                fac10(lay,icol) * absa(ind0+1,ig) + &
                fac01(lay,icol) * absa(ind1,  ig) + &
                fac11(lay,icol) * absa(ind1+1,ig))  &
               + tauself + taufor &
               + adjcolco2 * absco2 &
               + colcfc11(lay,icol) * cfc11adj(ig) &
               + colcfc12(lay,icol) * cfc12(ig)
            pfracs(lay,ngs5+ig,icol) = fracrefa(ig)
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         do ig = 1,ng6
            taug(lay,ngs5+ig,icol) = 0.0  &
               + colcfc11(lay,icol) * cfc11adj(ig) &
               + colcfc12(lay,icol) * cfc12(ig)
            pfracs(lay,ngs5+ig,icol) = fracrefa(ig)
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb6


   !-------------------------------------------
   subroutine taugb7 (ncol, nlay, taug, pfracs)
   !-------------------------------------------
   !
   !     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
   !                            (high key - o3; high minor - co2)

   ! Compute the optical depth by interpolating in ln(pressure), 
   ! temperature, and appropriate species. Below laytrop, the water
   ! vapor self-continuum and foreign continuum is interpolated 
   ! (in temperature) separately. 

      use parrrtm, only : ngs6
      use rrlw_ref, only : chi_mls
      use rrlw_kg07

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmco2, jpl
      real :: speccomb, specparm, specmult, fs
      real :: speccomb1, specparm1, specmult1, fs1
      real :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real :: p, p4, fk0, fk1, fk2
      real :: fac000, fac100, fac200, fac010, fac110, fac210
      real :: fac001, fac101, fac201, fac011, fac111, fac211
      real :: tauself, taufor, co2m1, co2m2, absco2
      real :: chi_co2, ratco2, adjfac, adjcolco2
      real :: refrat_planck_a, refrat_m_a
      real :: tau_major, tau_major1

      ! Minor gas mapping levels:
      !     lower - co2, p = 706.2620 mb, t = 278.94 K
      !     upper - co2, p = 12.9350 mb, t = 234.01 K

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower atmosphere.

      ! P = 706.2620 mb
      refrat_planck_a = chi_mls(1,3)/chi_mls(3,3)

      ! P = 706.2720 mb
      refrat_m_a = refrat_planck_a

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         speccomb = colh2o(lay,icol) + rat_h2oo3(lay,icol)*colo3(lay,icol)
         specparm = colh2o(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(lay,icol) + rat_h2oo3_1(lay,icol)*colo3(lay,icol)
         specparm1 = colh2o(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mco2 = colh2o(lay,icol) + refrat_m_a*colo3(lay,icol)
         specparm_mco2 = colh2o(lay,icol)/speccomb_mco2
         if (specparm_mco2 .ge. oneminus) specparm_mco2 = oneminus
         specmult_mco2 = 8. *specparm_mco2

         jmco2 = 1 + int(specmult_mco2)
         fmco2 = mod(specmult_mco2,1.0 )

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(lay,icol)/(coldry(lay,icol))
         ratco2 = 1.e20*chi_co2/chi_mls(2,jp(lay,icol)+1)
         if (ratco2 .gt. 3.0 ) then
            adjfac = 3.0 +(ratco2-3.0 )**0.79 
            adjcolco2 = adjfac*chi_mls(2,jp(lay,icol)+1)*coldry(lay,icol)*1.e-20 
         else
            adjcolco2 = colco2(lay,icol)
         endif

         speccomb_planck = colh2o(lay,icol)+refrat_planck_a*colo3(lay,icol)
         specparm_planck = colh2o(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(7) + js
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(7) + js1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else
            fac000 = (1.  - fs) * fac00(lay,icol)
            fac010 = (1.  - fs) * fac10(lay,icol)
            fac100 = fs * fac00(lay,icol)
            fac110 = fs * fac10(lay,icol)
         endif
         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else
            fac001 = (1.  - fs1) * fac01(lay,icol)
            fac011 = (1.  - fs1) * fac11(lay,icol)
            fac101 = fs1 * fac01(lay,icol)
            fac111 = fs1 * fac11(lay,icol)
         endif

         do ig = 1,ng7
            tauself = selffac(lay,icol)* (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
               (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
            co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
               (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(lay,icol) * (co2m2 - co2m1)

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac200 * absa(ind0+2, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig) + &
                   fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                  (fac200 * absa(ind0-1, ig) + &
                   fac100 * absa(ind0,   ig) + &
                   fac000 * absa(ind0+1, ig) + &
                   fac210 * absa(ind0+8, ig) + &
                   fac110 * absa(ind0+9, ig) + &
                   fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac201 * absa(ind1+2, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig) + &
                   fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                  (fac201 * absa(ind1-1, ig) + &
                   fac101 * absa(ind1,   ig) + &
                   fac001 * absa(ind1+1, ig) + &
                   fac211 * absa(ind1+8, ig) + &
                   fac111 * absa(ind1+9, ig) + &
                   fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) +  &
                   fac101 * absa(ind1+1, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs6+ig,icol) = &
               tau_major + tau_major1 + tauself + taufor &
                 + adjcolco2 * absco2
            pfracs(lay,ngs6+ig,icol) = fracrefa(ig,jpl) + fpl * &
               (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(lay,icol)/(coldry(lay,icol))
         ratco2 = 1.e20*chi_co2/chi_mls(2,jp(lay,icol)+1)
         if (ratco2 .gt. 3.0 ) then
            adjfac = 2.0 +(ratco2-2.0 )**0.79 
            adjcolco2 = adjfac*chi_mls(2,jp(lay,icol)+1)*coldry(lay,icol)*1.e-20 
         else
            adjcolco2 = colco2(lay,icol)
         endif

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(7) + 1
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(7) + 1
         indm = indminor(lay,icol)

         do ig = 1,ng7
            absco2 = kb_mco2(indm,ig) + minorfrac(lay,icol) * &
               (kb_mco2(indm+1,ig) - kb_mco2(indm,ig))
            taug(lay,ngs6+ig,icol) = colo3(lay,icol) * &
               (fac00(lay,icol) * absb(ind0,  ig) + &
                fac10(lay,icol) * absb(ind0+1,ig) + &
                fac01(lay,icol) * absb(ind1,  ig) + &
                fac11(lay,icol) * absb(ind1+1,ig))  &
               + adjcolco2 * absco2
            pfracs(lay,ngs6+ig,icol) = fracrefb(ig)
         enddo

! Empirical modification to code to improve stratospheric cooling rates
! for o3.  Revised to apply weighting for g-point reduction in this band.

         taug(lay,ngs6+6,icol)=taug(lay,ngs6+6,icol)*0.92 
         taug(lay,ngs6+7,icol)=taug(lay,ngs6+7,icol)*0.88 
         taug(lay,ngs6+8,icol)=taug(lay,ngs6+8,icol)*1.07 
         taug(lay,ngs6+9,icol)=taug(lay,ngs6+9,icol)*1.1 
         taug(lay,ngs6+10,icol)=taug(lay,ngs6+10,icol)*0.99 
         taug(lay,ngs6+11,icol)=taug(lay,ngs6+11,icol)*0.855 

      end do  ! upper layers

      end do  ! column

   end subroutine taugb7


   !-------------------------------------------
   subroutine taugb8 (ncol, nlay, taug, pfracs)
   !-------------------------------------------
   !
   !     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
   !                             (high key - o3; high minor - co2, n2o)

   ! Compute the optical depth by interpolating in ln(pressure) and 
   ! temperature, and appropriate species. Below laytrop, the water vapor 
   ! self-continuum and foreign continuum is interpolated (in temperature) 
   ! separately.

      use parrrtm, only : ngs7
      use rrlw_ref, only : chi_mls
      use rrlw_kg08

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, indm, ig
      real :: tauself, taufor, absco2, abso3, absn2o
      real :: chi_co2, ratco2, adjfac, adjcolco2

      ! Minor gas mapping levels:
      !     lower - co2, p = 1053.63 mb, t = 294.2 K
      !     lower - o3,  p = 317.348 mb, t = 240.77 K
      !     lower - n2o, p = 706.2720 mb, t = 278.94 K
      !     lower - cfc12,cfc11
      !     upper - co2, p = 35.1632 mb, t = 223.28 K
      !     upper - n2o, p = 8.716e-2 mb, t = 226.03 K

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(lay,icol)/(coldry(lay,icol))
         ratco2 = 1.e20 *chi_co2/chi_mls(2,jp(lay,icol)+1)
         if (ratco2 .gt. 3.0 ) then
            adjfac = 2.0 +(ratco2-2.0 )**0.65 
            adjcolco2 = adjfac*chi_mls(2,jp(lay,icol)+1)*coldry(lay,icol)*1.e-20 
         else
            adjcolco2 = colco2(lay,icol)
         endif

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(8) + 1
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(8) + 1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)

         do ig = 1,ng8
            tauself = selffac(lay,icol) * (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig)))
            absco2 = (ka_mco2(indm,ig) + minorfrac(lay,icol) * &
               (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            abso3 = (ka_mo3(indm,ig) + minorfrac(lay,icol) * &
               (ka_mo3(indm+1,ig) - ka_mo3(indm,ig)))
            absn2o = (ka_mn2o(indm,ig) + minorfrac(lay,icol) * &
               (ka_mn2o(indm+1,ig) - ka_mn2o(indm,ig)))
            taug(lay,ngs7+ig,icol) = colh2o(lay,icol) * &
               (fac00(lay,icol) * absa(ind0,ig)   + &
                fac10(lay,icol) * absa(ind0+1,ig) + &
                fac01(lay,icol) * absa(ind1,ig)   + &
                fac11(lay,icol) * absa(ind1+1,ig))  &
               + tauself + taufor &
               + adjcolco2 * absco2 &
               + colo3(lay,icol) * abso3 &
               + coln2o(lay,icol) * absn2o &
               + colcfc12(lay,icol) * cfc12(ig) &
               + colcfc22(lay,icol) * cfc22adj(ig)
            pfracs(lay,ngs7+ig,icol) = fracrefa(ig)
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(lay,icol)/coldry(lay,icol)
         ratco2 = 1.e20 *chi_co2/chi_mls(2,jp(lay,icol)+1)
         if (ratco2 .gt. 3.0 ) then
            adjfac = 2.0 +(ratco2-2.0 )**0.65 
            adjcolco2 = adjfac*chi_mls(2,jp(lay,icol)+1) * coldry(lay,icol)*1.e-20 
         else
            adjcolco2 = colco2(lay,icol)
         endif

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(8) + 1
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(8) + 1
         indm = indminor(lay,icol)

         do ig = 1,ng8
            absco2 = (kb_mco2(indm,ig) + minorfrac(lay,icol) * &
               (kb_mco2(indm+1,ig) - kb_mco2(indm,ig)))
            absn2o = (kb_mn2o(indm,ig) + minorfrac(lay,icol) * &
               (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig)))
            taug(lay,ngs7+ig,icol) = colo3(lay,icol) * &
               (fac00(lay,icol) * absb(ind0,  ig) + &
                fac10(lay,icol) * absb(ind0+1,ig) + &
                fac01(lay,icol) * absb(ind1,  ig) + &
                fac11(lay,icol) * absb(ind1+1,ig)) &
               + adjcolco2 * absco2 &
               + coln2o(lay,icol) * absn2o & 
               + colcfc12(lay,icol) * cfc12(ig) &
               + colcfc22(lay,icol) * cfc22adj(ig)
            pfracs(lay,ngs7+ig,icol) = fracrefb(ig)
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb8


   !-------------------------------------------
   subroutine taugb9 (ncol, nlay, taug, pfracs)
   !-------------------------------------------
   !
   !     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
   !                             (high key - ch4; high minor - n2o)

   ! Compute the optical depth by interpolating in ln(pressure), 
   ! temperature, and appropriate species. Below laytrop, the water
   ! vapor self-continuum and foreign continuum is interpolated 
   ! (in temperature) separately.  

      use parrrtm, only : ngs8
      use rrlw_ref, only : chi_mls
      use rrlw_kg09

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer  :: icol, lay, ind0, ind1, inds, indf, indm, ig
      integer  :: js, js1, jmn2o, jpl
      real :: speccomb, specparm, specmult, fs
      real :: speccomb1, specparm1, specmult1, fs1
      real :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o
      real :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real :: p, p4, fk0, fk1, fk2
      real :: fac000, fac100, fac200, fac010, fac110, fac210
      real :: fac001, fac101, fac201, fac011, fac111, fac211
      real :: tauself, taufor, n2om1, n2om2, absn2o
      real :: chi_n2o, ratn2o, adjfac, adjcoln2o
      real :: refrat_planck_a, refrat_m_a
      real :: tau_major, tau_major1

      ! Minor gas mapping levels:
      !     lower - n2o, p = 706.272 mb, t = 278.94 K
      !     upper - n2o, p = 95.58 mb, t = 215.7 K

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower atmosphere.

      ! P = 212 mb
      refrat_planck_a = chi_mls(1,9)/chi_mls(6,9)

      ! P = 706.272 mb 
      refrat_m_a = chi_mls(1,3)/chi_mls(6,3)

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         speccomb = colh2o(lay,icol) + rat_h2och4(lay,icol)*colch4(lay,icol)
         specparm = colh2o(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(lay,icol) + rat_h2och4_1(lay,icol)*colch4(lay,icol)
         specparm1 = colh2o(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mn2o = colh2o(lay,icol) + refrat_m_a*colch4(lay,icol)
         specparm_mn2o = colh2o(lay,icol)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminus) specparm_mn2o = oneminus
         specmult_mn2o = 8. *specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0 )

!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(lay,icol)/(coldry(lay,icol))
         ratn2o = 1.e20 *chi_n2o/chi_mls(4,jp(lay,icol)+1)
         if (ratn2o .gt. 1.5 ) then
            adjfac = 0.5 +(ratn2o-0.5 )**0.65 
            adjcoln2o = adjfac*chi_mls(4,jp(lay,icol)+1)*coldry(lay,icol)*1.e-20 
         else
            adjcoln2o = coln2o(lay,icol)
         endif

         speccomb_planck = colh2o(lay,icol)+refrat_planck_a*colch4(lay,icol)
         specparm_planck = colh2o(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(9) + js
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(9) + js1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else
            fac000 = (1.  - fs) * fac00(lay,icol)
            fac010 = (1.  - fs) * fac10(lay,icol)
            fac100 = fs * fac00(lay,icol)
            fac110 = fs * fac10(lay,icol)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else
            fac001 = (1.  - fs1) * fac01(lay,icol)
            fac011 = (1.  - fs1) * fac11(lay,icol)
            fac101 = fs1 * fac01(lay,icol)
            fac111 = fs1 * fac11(lay,icol)
         endif

         do ig = 1,ng9
            tauself = selffac(lay,icol)* (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
               (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
            n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
               (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(lay,icol) * (n2om2 - n2om1)

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac200 * absa(ind0+2, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig) + &
                   fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                  (fac200 * absa(ind0-1, ig) + &
                   fac100 * absa(ind0,   ig) + &
                   fac000 * absa(ind0+1, ig) + &
                   fac210 * absa(ind0+8, ig) + &
                   fac110 * absa(ind0+9, ig) + &
                   fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + & 
                   fac101 * absa(ind1+1, ig) + &
                   fac201 * absa(ind1+2, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig) + &
                   fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                  (fac201 * absa(ind1-1, ig) + &
                   fac101 * absa(ind1,   ig) + &
                   fac001 * absa(ind1+1, ig) + &
                   fac211 * absa(ind1+8, ig) + &
                   fac111 * absa(ind1+9, ig) + &
                   fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs8+ig,icol) = &
               tau_major + tau_major1 + tauself + taufor &
                 + adjcoln2o * absn2o
            pfracs(lay,ngs8+ig,icol) = fracrefa(ig,jpl) + fpl * &
               (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(lay,icol)/(coldry(lay,icol))
         ratn2o = 1.e20 *chi_n2o/chi_mls(4,jp(lay,icol)+1)
         if (ratn2o .gt. 1.5 ) then
            adjfac = 0.5 +(ratn2o-0.5 )**0.65 
            adjcoln2o = adjfac*chi_mls(4,jp(lay,icol)+1)*coldry(lay,icol)*1.e-20 
         else
            adjcoln2o = coln2o(lay,icol)
         endif

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(9) + 1
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(9) + 1
         indm = indminor(lay,icol)

         do ig = 1,ng9
            absn2o = kb_mn2o(indm,ig) + minorfrac(lay,icol) * &
               (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig))
            taug(lay,ngs8+ig,icol) = colch4(lay,icol) * &
               (fac00(lay,icol) * absb(ind0,  ig) + &
                fac10(lay,icol) * absb(ind0+1,ig) + &
                fac01(lay,icol) * absb(ind1,  ig) + &
                fac11(lay,icol) * absb(ind1+1,ig)) &
               + adjcoln2o * absn2o
            pfracs(lay,ngs8+ig,icol) = fracrefb(ig)
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb9


   !--------------------------------------------
   subroutine taugb10 (ncol, nlay, taug, pfracs)
   !--------------------------------------------
   !
   !     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)

   ! Compute the optical depth by interpolating in ln(pressure) and 
   ! temperature. Below laytrop, the water vapor self-continuum and
   ! foreign continuum is interpolated (in temperature) separately.

      use parrrtm, only : ngs9
      use rrlw_kg10

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, ig
      real :: tauself, taufor

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(10) + 1
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(10) + 1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)

         do ig = 1,ng10
            tauself = selffac(lay,icol) * (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            taug(lay,ngs9+ig,icol) = colh2o(lay,icol) * &
               (fac00(lay,icol) * absa(ind0,  ig) + &
                fac10(lay,icol) * absa(ind0+1,ig) + &
                fac01(lay,icol) * absa(ind1,  ig) + &
                fac11(lay,icol) * absa(ind1+1,ig))  &
               + tauself + taufor
            pfracs(lay,ngs9+ig,icol) = fracrefa(ig)
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(10) + 1
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(10) + 1
         indf = indfor(lay,icol)

         do ig = 1,ng10
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            taug(lay,ngs9+ig,icol) = colh2o(lay,icol) * &
               (fac00(lay,icol) * absb(ind0,  ig) + &
                fac10(lay,icol) * absb(ind0+1,ig) + &
                fac01(lay,icol) * absb(ind1,  ig) + &
                fac11(lay,icol) * absb(ind1+1,ig))  &
               + taufor
            pfracs(lay,ngs9+ig,icol) = fracrefb(ig)
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb10


   !--------------------------------------------
   subroutine taugb11 (ncol, nlay, taug, pfracs)
   !--------------------------------------------
   !
   !     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
   !                              (high key - h2o; high minor - o2)

   ! Compute the optical depth by interpolating in ln(pressure) and 
   ! temperature. Below laytrop, the water vapor self-continuum and
   ! foreign continuum is interpolated (in temperature) separately.

      use parrrtm, only : ngs10
      use rrlw_kg11

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, indm, ig
      real :: scaleo2, tauself, taufor, tauo2

      ! Minor gas mapping levels:
      !     lower - o2, p = 706.2720 mb, t = 278.94 K
      !     upper - o2, p = 4.758820 mb, t = 250.85 K

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(11) + 1
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(11) + 1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)
         scaleo2 = colo2(lay,icol)*scaleminor(lay,icol)
         do ig = 1,ng11
            tauself = selffac(lay,icol) * (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig)))
            tauo2 = scaleo2 * (ka_mo2(indm,ig) + minorfrac(lay,icol) * &
               (ka_mo2(indm+1,ig) - ka_mo2(indm,ig)))
            taug(lay,ngs10+ig,icol) = colh2o(lay,icol) * &
               (fac00(lay,icol) * absa(ind0,  ig) + &
                fac10(lay,icol) * absa(ind0+1,ig) + &
                fac01(lay,icol) * absa(ind1,  ig) + &
                fac11(lay,icol) * absa(ind1+1,ig))  &
               + tauself + taufor + tauo2
            pfracs(lay,ngs10+ig,icol) = fracrefa(ig)
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(11) + 1
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(11) + 1
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)
         scaleo2 = colo2(lay,icol)*scaleminor(lay,icol)
         do ig = 1,ng11
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            tauo2 = scaleo2 * (kb_mo2(indm,ig) + minorfrac(lay,icol) * &
               (kb_mo2(indm+1,ig) - kb_mo2(indm,ig)))
            taug(lay,ngs10+ig,icol) = colh2o(lay,icol) * &
               (fac00(lay,icol) * absb(ind0,  ig) + &
                fac10(lay,icol) * absb(ind0+1,ig) + &
                fac01(lay,icol) * absb(ind1,  ig) + &
                fac11(lay,icol) * absb(ind1+1,ig))  &
               + taufor + tauo2
            pfracs(lay,ngs10+ig,icol) = fracrefb(ig)
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb11


   !--------------------------------------------
   subroutine taugb12 (ncol, nlay, taug, pfracs)
   !--------------------------------------------
   !
   !     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)

   ! Compute the optical depth by interpolating in ln(pressure), 
   ! temperature, and appropriate species. Below laytrop, the water
   ! vapor self-continuum adn foreign continuum is interpolated 
   ! (in temperature) separately.  

      use parrrtm, only : ngs11
      use rrlw_ref, only : chi_mls
      use rrlw_kg12

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, ig
      integer :: js, js1, jpl
      real :: speccomb, specparm, specmult, fs
      real :: speccomb1, specparm1, specmult1, fs1
      real :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real :: p, p4, fk0, fk1, fk2
      real :: fac000, fac100, fac200, fac010, fac110, fac210
      real :: fac001, fac101, fac201, fac011, fac111, fac211
      real :: tauself, taufor
      real :: refrat_planck_a
      real :: tau_major, tau_major1

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower atmosphere.

      ! P = 174.164 mb 
      refrat_planck_a = chi_mls(1,10)/chi_mls(2,10)

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         speccomb = colh2o(lay,icol) + rat_h2oco2(lay,icol)*colco2(lay,icol)
         specparm = colh2o(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(lay,icol) + rat_h2oco2_1(lay,icol)*colco2(lay,icol)
         specparm1 = colh2o(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_planck = colh2o(lay,icol)+refrat_planck_a*colco2(lay,icol)
         specparm_planck = colh2o(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(12) + js
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(12) + js1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else
            fac000 = (1.  - fs) * fac00(lay,icol)
            fac010 = (1.  - fs) * fac10(lay,icol)
            fac100 = fs * fac00(lay,icol)
            fac110 = fs * fac10(lay,icol)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else
            fac001 = (1.  - fs1) * fac01(lay,icol)
            fac011 = (1.  - fs1) * fac11(lay,icol)
            fac101 = fs1 * fac01(lay,icol)
            fac111 = fs1 * fac11(lay,icol)
         endif

         do ig = 1,ng12
            tauself = selffac(lay,icol) * (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac200 * absa(ind0+2, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig) + &
                   fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                  (fac200 * absa(ind0-1, ig) + &
                   fac100 * absa(ind0,   ig) + &
                   fac000 * absa(ind0+1, ig) + &
                   fac210 * absa(ind0+8, ig) + &
                   fac110 * absa(ind0+9, ig) + &
                   fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac201 * absa(ind1+2, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig) + &
                   fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                  (fac201 * absa(ind1-1, ig) + &
                   fac101 * absa(ind1,   ig) + &
                   fac001 * absa(ind1+1, ig) + &
                   fac211 * absa(ind1+8, ig) + &
                   fac111 * absa(ind1+9, ig) + &
                   fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs11+ig,icol) = &
               tau_major + tau_major1 + tauself + taufor
            pfracs(lay,ngs11+ig,icol) = fracrefa(ig,jpl) + fpl * &
               (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
   
      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         do ig = 1,ng12
            taug(lay,ngs11+ig,icol) = 0.0 
            pfracs(lay,ngs11+ig,icol) = 0.0 
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb12


   !--------------------------------------------
   subroutine taugb13 (ncol, nlay, taug, pfracs)
   !--------------------------------------------
   !
   !     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)

   ! Compute the optical depth by interpolating in ln(pressure), 
   ! temperature, and appropriate species. Below laytrop, the water
   ! vapor self-continuum and foreign continuum is interpolated 
   ! (in temperature) separately.  

      use parrrtm, only : ngs12
      use rrlw_ref, only : chi_mls
      use rrlw_kg13

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmco2, jmco, jpl
      real :: speccomb, specparm, specmult, fs
      real :: speccomb1, specparm1, specmult1, fs1
      real :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real :: speccomb_mco, specparm_mco, specmult_mco, fmco
      real :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real :: p, p4, fk0, fk1, fk2
      real :: fac000, fac100, fac200, fac010, fac110, fac210
      real :: fac001, fac101, fac201, fac011, fac111, fac211
      real :: tauself, taufor, co2m1, co2m2, absco2 
      real :: com1, com2, absco, abso3
      real :: chi_co2, ratco2, adjfac, adjcolco2
      real :: refrat_planck_a, refrat_m_a, refrat_m_a3
      real :: tau_major, tau_major1

      ! Minor gas mapping levels :
      !     lower - co2, p = 1053.63 mb, t = 294.2 K
      !     lower - co, p = 706 mb, t = 278.94 K
      !     upper - o3, p = 95.5835 mb, t = 215.7 K

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower atmosphere.

      ! P = 473.420 mb (Level 5)
      refrat_planck_a = chi_mls(1,5)/chi_mls(4,5)

      ! P = 1053. (Level 1)
      refrat_m_a = chi_mls(1,1)/chi_mls(4,1)

      ! P = 706. (Level 3)
      refrat_m_a3 = chi_mls(1,3)/chi_mls(4,3)

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         speccomb = colh2o(lay,icol) + rat_h2on2o(lay,icol)*coln2o(lay,icol)
         specparm = colh2o(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(lay,icol) + rat_h2on2o_1(lay,icol)*coln2o(lay,icol)
         specparm1 = colh2o(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mco2 = colh2o(lay,icol) + refrat_m_a*coln2o(lay,icol)
         specparm_mco2 = colh2o(lay,icol)/speccomb_mco2
         if (specparm_mco2 .ge. oneminus) specparm_mco2 = oneminus
         specmult_mco2 = 8. *specparm_mco2
         jmco2 = 1 + int(specmult_mco2)
         fmco2 = mod(specmult_mco2,1.0 )

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(lay,icol)/(coldry(lay,icol))
         ratco2 = 1.e20 *chi_co2/3.55e-4 
         if (ratco2 .gt. 3.0 ) then
            adjfac = 2.0 +(ratco2-2.0 )**0.68 
            adjcolco2 = adjfac*3.55e-4*coldry(lay,icol)*1.e-20 
         else
            adjcolco2 = colco2(lay,icol)
         endif

         speccomb_mco = colh2o(lay,icol) + refrat_m_a3*coln2o(lay,icol)
         specparm_mco = colh2o(lay,icol)/speccomb_mco
         if (specparm_mco .ge. oneminus) specparm_mco = oneminus
         specmult_mco = 8. *specparm_mco
         jmco = 1 + int(specmult_mco)
         fmco = mod(specmult_mco,1.0 )

         speccomb_planck = colh2o(lay,icol)+refrat_planck_a*coln2o(lay,icol)
         specparm_planck = colh2o(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(13) + js
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(13) + js1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else
            fac000 = (1.  - fs) * fac00(lay,icol)
            fac010 = (1.  - fs) * fac10(lay,icol)
            fac100 = fs * fac00(lay,icol)
            fac110 = fs * fac10(lay,icol)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else
            fac001 = (1.  - fs1) * fac01(lay,icol)
            fac011 = (1.  - fs1) * fac11(lay,icol)
            fac101 = fs1 * fac01(lay,icol)
            fac111 = fs1 * fac11(lay,icol)
         endif

         do ig = 1,ng13
            tauself = selffac(lay,icol) * (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
               (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
            co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
               (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(lay,icol) * (co2m2 - co2m1)
            com1 = ka_mco(jmco,indm,ig) + fmco * &
               (ka_mco(jmco+1,indm,ig) - ka_mco(jmco,indm,ig))
            com2 = ka_mco(jmco,indm+1,ig) + fmco * &
               (ka_mco(jmco+1,indm+1,ig) - ka_mco(jmco,indm+1,ig))
            absco = com1 + minorfrac(lay,icol) * (com2 - com1)

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac200 * absa(ind0+2, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig) + &
                   fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                  (fac200 * absa(ind0-1, ig) + &
                   fac100 * absa(ind0,   ig) + &
                   fac000 * absa(ind0+1, ig) + &
                   fac210 * absa(ind0+8, ig) + &
                   fac110 * absa(ind0+9, ig) + &
                   fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac201 * absa(ind1+2, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig) + &
                   fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                  (fac201 * absa(ind1-1, ig) + &
                   fac101 * absa(ind1,   ig) + &
                   fac001 * absa(ind1+1, ig) + &
                   fac211 * absa(ind1+8, ig) + &
                   fac111 * absa(ind1+9, ig) + &
                   fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs12+ig,icol) = &
               tau_major + tau_major1 + tauself + taufor &
                 + adjcolco2 * absco2 &
                 + colco(lay,icol) * absco
            pfracs(lay,ngs12+ig,icol) = fracrefa(ig,jpl) + fpl * &
               (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         indm = indminor(lay,icol)
         do ig = 1,ng13
            abso3 = kb_mo3(indm,ig) + minorfrac(lay,icol) * &
               (kb_mo3(indm+1,ig) - kb_mo3(indm,ig))
            taug(lay,ngs12+ig,icol) = colo3(lay,icol) * abso3
            pfracs(lay,ngs12+ig,icol) = fracrefb(ig)
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb13


   !---------------------------------------------
   subroutine taugb14 (ncol, nlay , taug, pfracs)
   !---------------------------------------------
   !
   !     band 14:  2250-2380 cm-1 (low - co2; high - co2)

   ! Compute the optical depth by interpolating in ln(pressure) and 
   ! temperature. Below laytrop, the water vapor self-continuum 
   ! and foreign continuum is interpolated (in temperature) separately.  

      use parrrtm, only : ngs13
      use rrlw_kg14

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, ig
      real :: tauself, taufor

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(14) + 1
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(14) + 1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         do ig = 1,ng14
            tauself = selffac(lay,icol) * (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            taug(lay,ngs13+ig,icol) = colco2(lay,icol) * &
               (fac00(lay,icol) * absa(ind0,  ig) + &
                fac10(lay,icol) * absa(ind0+1,ig) + &
                fac01(lay,icol) * absa(ind1,  ig) + &
                fac11(lay,icol) * absa(ind1+1,ig))  &
               + tauself + taufor
            pfracs(lay,ngs13+ig,icol) = fracrefa(ig)
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(14) + 1
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(14) + 1
         do ig = 1,ng14
            taug(lay,ngs13+ig,icol) = colco2(lay,icol) * &
               (fac00(lay,icol) * absb(ind0,  ig) + &
                fac10(lay,icol) * absb(ind0+1,ig) + &
                fac01(lay,icol) * absb(ind1,  ig) + &
                fac11(lay,icol) * absb(ind1+1,ig))
            pfracs(lay,ngs13+ig,icol) = fracrefb(ig)
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb14


   !--------------------------------------------
   subroutine taugb15 (ncol, nlay, taug, pfracs)
   !--------------------------------------------
   !
   !     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
   !                              (high - nothing)

   ! Compute the optical depth by interpolating in ln(pressure), 
   ! temperature, and appropriate species. Below laytrop, the water
   ! vapor self-continuum and foreign continuum is interpolated 
   ! (in temperature) separately.  

      use parrrtm, only : ngs14
      use rrlw_ref, only : chi_mls
      use rrlw_kg15

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmn2, jpl
      real :: speccomb, specparm, specmult, fs
      real :: speccomb1, specparm1, specmult1, fs1
      real :: speccomb_mn2, specparm_mn2, specmult_mn2, fmn2
      real :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real :: p, p4, fk0, fk1, fk2
      real :: fac000, fac100, fac200, fac010, fac110, fac210
      real :: fac001, fac101, fac201, fac011, fac111, fac211
      real :: scalen2, tauself, taufor, n2m1, n2m2, taun2 
      real :: refrat_planck_a, refrat_m_a
      real :: tau_major, tau_major1

      ! Minor gas mapping levels: 
      !     Lower - Nitrogen Continuum, P = 1053., T = 294.

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower atmosphere.

      ! P = 1053. mb (Level 1)
      refrat_planck_a = chi_mls(4,1)/chi_mls(2,1)
      refrat_m_a = refrat_planck_a

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         speccomb = coln2o(lay,icol) + rat_n2oco2(lay,icol)*colco2(lay,icol)
         specparm = coln2o(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = coln2o(lay,icol) + rat_n2oco2_1(lay,icol)*colco2(lay,icol)
         specparm1 = coln2o(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mn2 = coln2o(lay,icol) + refrat_m_a*colco2(lay,icol)
         specparm_mn2 = coln2o(lay,icol)/speccomb_mn2
         if (specparm_mn2 .ge. oneminus) specparm_mn2 = oneminus
         specmult_mn2 = 8. *specparm_mn2
         jmn2 = 1 + int(specmult_mn2)
         fmn2 = mod(specmult_mn2,1.0 )

         speccomb_planck = coln2o(lay,icol)+refrat_planck_a*colco2(lay,icol)
         specparm_planck = coln2o(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(15) + js
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(15) + js1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)
         indm = indminor(lay,icol)
         
         scalen2 = colbrd(lay,icol)*scaleminor(lay,icol)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else
            fac000 = (1.  - fs) * fac00(lay,icol)
            fac010 = (1.  - fs) * fac10(lay,icol)
            fac100 = fs * fac00(lay,icol)
            fac110 = fs * fac10(lay,icol)
         endif
         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else
            fac001 = (1.  - fs1) * fac01(lay,icol)
            fac011 = (1.  - fs1) * fac11(lay,icol)
            fac101 = fs1 * fac01(lay,icol)
            fac111 = fs1 * fac11(lay,icol)
         endif

         do ig = 1,ng15
            tauself = selffac(lay,icol)* (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 
            n2m1 = ka_mn2(jmn2,indm,ig) + fmn2 * &
               (ka_mn2(jmn2+1,indm,ig) - ka_mn2(jmn2,indm,ig))
            n2m2 = ka_mn2(jmn2,indm+1,ig) + fmn2 * &
               (ka_mn2(jmn2+1,indm+1,ig) - ka_mn2(jmn2,indm+1,ig))
            taun2 = scalen2 * (n2m1 + minorfrac(lay,icol) * (n2m2 - n2m1))

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac200 * absa(ind0+2, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig) + &
                   fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                  (fac200 * absa(ind0-1, ig) + &
                   fac100 * absa(ind0,   ig) + &
                   fac000 * absa(ind0+1, ig) + &
                   fac210 * absa(ind0+8, ig) + &
                   fac110 * absa(ind0+9, ig) + &
                   fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig))
            endif 

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac201 * absa(ind1+2, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig) + &
                   fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                  (fac201 * absa(ind1-1, ig) + &
                   fac101 * absa(ind1,   ig) + &
                   fac001 * absa(ind1+1, ig) + &
                   fac211 * absa(ind1+8, ig) + &
                   fac111 * absa(ind1+9, ig) + &
                   fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs14+ig,icol) = &
               tau_major + tau_major1 + tauself + taufor + taun2
            pfracs(lay,ngs14+ig,icol) = fracrefa(ig,jpl) + fpl * &
               (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
    
      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         do ig = 1,ng15
            taug(lay,ngs14+ig,icol) = 0.0 
            pfracs(lay,ngs14+ig,icol) = 0.0 
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb15


   !--------------------------------------------
   subroutine taugb16 (ncol, nlay, taug, pfracs)
   !--------------------------------------------
   !
   !     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)

   ! Compute the optical depth by interpolating in ln(pressure), 
   ! temperature,and appropriate species. Below laytrop, the water
   ! vapor self-continuum and foreign continuum is interpolated 
   ! (in temperature) separately.  

      use parrrtm, only : ngs15
      use rrlw_ref, only : chi_mls
      use rrlw_kg16

      integer, intent(in)    :: ncol, nlay
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)
      real,    intent(inout) :: pfracs(nlay,ngptlw,ncol)

      integer :: icol, lay, ind0, ind1, inds, indf, ig
      integer :: js, js1, jpl
      real :: speccomb, specparm, specmult, fs
      real :: speccomb1, specparm1, specmult1, fs1
      real :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real :: p, p4, fk0, fk1, fk2
      real :: fac000, fac100, fac200, fac010, fac110, fac210
      real :: fac001, fac101, fac201, fac011, fac111, fac211
      real :: tauself, taufor
      real :: refrat_planck_a
      real :: tau_major, tau_major1

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower atmosphere.

      ! P = 387. mb (Level 6)
      refrat_planck_a = chi_mls(1,6)/chi_mls(6,6)

      do icol = 1,ncol

      ! lower atmosphere
      do lay = 1,laytrop(icol)

         speccomb = colh2o(lay,icol) + rat_h2och4(lay,icol)*colch4(lay,icol)
         specparm = colh2o(lay,icol)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(lay,icol) + rat_h2och4_1(lay,icol)*colch4(lay,icol)
         specparm1 = colh2o(lay,icol)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_planck = colh2o(lay,icol)+refrat_planck_a*colch4(lay,icol)
         specparm_planck = colh2o(lay,icol)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(lay,icol)-1)*5+(jt(lay,icol)-1))*nspa(16) + js
         ind1 = (jp(lay,icol)*5+(jt1(lay,icol)-1))*nspa(16) + js1
         inds = indself(lay,icol)
         indf = indfor(lay,icol)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay,icol)
            fac100 = fk1*fac00(lay,icol)
            fac200 = fk2*fac00(lay,icol)
            fac010 = fk0*fac10(lay,icol)
            fac110 = fk1*fac10(lay,icol)
            fac210 = fk2*fac10(lay,icol)
         else
            fac000 = (1.  - fs) * fac00(lay,icol)
            fac010 = (1.  - fs) * fac10(lay,icol)
            fac100 = fs * fac00(lay,icol)
            fac110 = fs * fac10(lay,icol)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay,icol)
            fac101 = fk1*fac01(lay,icol)
            fac201 = fk2*fac01(lay,icol)
            fac011 = fk0*fac11(lay,icol)
            fac111 = fk1*fac11(lay,icol)
            fac211 = fk2*fac11(lay,icol)
         else
            fac001 = (1. - fs1) * fac01(lay,icol)
            fac011 = (1. - fs1) * fac11(lay,icol)
            fac101 = fs1 * fac01(lay,icol)
            fac111 = fs1 * fac11(lay,icol)
         endif

         do ig = 1,ng16
            tauself = selffac(lay,icol)* (selfref(inds,ig) + selffrac(lay,icol) * &
               (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay,icol) * (forref(indf,ig) + forfrac(lay,icol) * &
               (forref(indf+1,ig) - forref(indf,ig))) 

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac200 * absa(ind0+2, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig) + &
                   fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                  (fac200 * absa(ind0-1, ig) + &
                   fac100 * absa(ind0,   ig) + &
                   fac000 * absa(ind0+1, ig) + &
                   fac210 * absa(ind0+8, ig) + &
                   fac110 * absa(ind0+9, ig) + &
                   fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125) then
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac201 * absa(ind1+2, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig) + &
                   fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875) then
               tau_major1 = speccomb1 * &
                  (fac201 * absa(ind1-1, ig) + &
                   fac101 * absa(ind1,   ig) + &
                   fac001 * absa(ind1+1, ig) + &
                   fac211 * absa(ind1+8, ig) + &
                   fac111 * absa(ind1+9, ig) + &
                   fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                  (fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs15+ig,icol) = &
               tau_major + tau_major1 + tauself + taufor
            pfracs(lay,ngs15+ig,icol) = fracrefa(ig,jpl) + fpl * &
               (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo

      end do  ! lower layers

      ! upper atmosphere
      do lay = laytrop(icol)+1, nlay

         ind0 = ((jp(lay,icol)-13)*5+(jt(lay,icol)-1))*nspb(16) + 1
         ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(16) + 1
         do ig = 1,ng16
            taug(lay,ngs15+ig,icol) = colch4(lay,icol) * &
               (fac00(lay,icol) * absb(ind0,  ig) + &
                fac10(lay,icol) * absb(ind0+1,ig) + &
                fac01(lay,icol) * absb(ind1,  ig) + &
                fac11(lay,icol) * absb(ind1+1,ig))
            pfracs(lay,ngs15+ig,icol) = fracrefb(ig)
         enddo

      end do  ! upper layers

      end do  ! column

   end subroutine taugb16


   ! ---------------------------------------------
   subroutine addAerosols (ncol, nlay, taua, taug)
   ! ---------------------------------------------

      integer, intent(in)    :: ncol, nlay
      real,    intent(in)    :: taua(nlay,nbndlw,ncol)
      real,    intent(inout) :: taug(nlay,ngptlw,ncol)

      integer :: icol, lay, ig, ibnd
     
      do icol = 1,ncol
         do ig = 1,ngptlw
            ibnd = ngb(ig)
            taug(:,ig,icol) = taug(:,ig,icol) + taua(:,ibnd,icol)
         end do
      end do

   end subroutine addAerosols

end module rrtmg_lw_taumol

