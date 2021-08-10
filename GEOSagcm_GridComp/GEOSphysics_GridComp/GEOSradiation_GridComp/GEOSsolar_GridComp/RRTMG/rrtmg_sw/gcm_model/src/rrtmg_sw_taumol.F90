!pmn use macro for linear interp?
!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! space saving function-like macros for linear interpolation
! of a 2D-varaible in first and second arguments as follows ...
#define LIN2_ARG1(VAR,I,J,FINT) (VAR(I,J) + FINT * (VAR(I+1,J)-VAR(I,J)))
#define LIN2_ARG2(VAR,I,J,FINT) (VAR(I,J) + FINT * (VAR(I,J+1)-VAR(I,J)))

module rrtmg_sw_taumol

   use parrrsw, only: ngptsw, jpband
   use rrsw_con, only: oneminus
   use rrsw_wvn, only: nspa, nspb

   implicit none

contains

   !----------------------------------------------------------------------------
   subroutine taumol_sw(pncol, ncol, nlayers, &
                        colh2o, colco2, colch4, colo2, colo3, colmol, &
                        laytrop, jp, jt, jt1, &
                        fac00, fac01, fac10, fac11, &
                        selffac, selffrac, indself, forfac, forfrac, indfor, &
                        isolvar, svar_f, svar_s, svar_i, &
                        svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                        ssi, sfluxzen, taug, taur)

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)       ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! Calculate gaseous optical depth and planck fractions for each spectral band.

      call taumol16(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol17(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol18(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol19(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol20(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol21(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol22(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol23(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol24(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol25(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol26(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol27(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol28(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

      call taumol29(pncol, ncol, nlayers, &
                    colh2o, colco2, colch4, colo2, colo3, colmol, &
                    laytrop, jp, jt, jt1, &
                    fac00, fac01, fac10, fac11, &
                    isolvar, svar_f, svar_s, svar_i, &
                    svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                    selffac, selffrac, indself, forfac, forfrac, indfor, &
                    ssi, sfluxzen, taug, taur)

   end subroutine

   !----------------------------------------------------------------------------
   subroutine taumol16(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng16
      use rrsw_kg16, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray, strrat1

      strrat1 = 252.131 
      layreffr = 18

      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  

      do icol = 1,ncol

         ! Lower atmosphere loop
         do lay = 1,nlayers
            if (lay <= laytrop(icol)) then 
               speccomb = colh2o(lay,icol) + strrat1*colch4(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               fac000 = (1. - fs) * fac00(lay,icol) 
               fac010 = (1. - fs) * fac10(lay,icol) 
               fac100 =       fs  * fac00(lay,icol) 
               fac110 =       fs  * fac10(lay,icol) 
               fac001 = (1. - fs) * fac01(lay,icol) 
               fac011 = (1. - fs) * fac11(lay,icol) 
               fac101 =       fs  * fac01(lay,icol) 
               fac111 =       fs  * fac11(lay,icol) 
               ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(16) + js
               ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(16) + js
               inds = indself(lay,icol) 
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng16
                  taug(lay,ig,icol) = speccomb * &
                     (fac000 * absa(ind0   ,ig) + &
                      fac100 * absa(ind0 +1,ig) + &
                      fac010 * absa(ind0 +9,ig) + &
                      fac110 * absa(ind0+10,ig) + &
                      fac001 * absa(ind1   ,ig) + &
                      fac101 * absa(ind1 +1,ig) + &
                      fac011 * absa(ind1 +9,ig) + &
                      fac111 * absa(ind1+10,ig)) + &
                     colh2o(lay,icol) * &
                     (selffac(lay,icol) * LIN2_ARG1(selfref,inds,ig,selffrac(lay,icol)) + &
                      forfac (lay,icol) * LIN2_ARG1( forref,indf,ig, forfrac(lay,icol)))
!                 ssa(lay,ig) = tauray / taug(lay,ig)
                  taur(lay,ig,icol) = tauray
    
               enddo
            end if
         enddo
      end do

      ! Upper atmosphere loop
      do icol = 1,ncol
         laysolfr = nlayers
         do lay = 1,nlayers
            if (lay > laytrop(icol)) then
            !do lay = laytrop(icol) +1, nlayers
               if (jp(lay-1,icol) < layreffr .and. jp(lay,icol) >= layreffr) then
                  laysolfr = lay
               end if
               ind0 = ((jp(lay,icol)-13)*5+(jt (lay,icol)-1))*nspb(16) + 1
               ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(16) + 1
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng16
                  taug(lay,ig,icol) = colch4(lay,icol) * &
                     (fac00(lay,icol) * absb(ind0  ,ig) + &
                      fac10(lay,icol) * absb(ind0+1,ig) + &
                      fac01(lay,icol) * absb(ind1  ,ig) + &
                      fac11(lay,icol) * absb(ind1+1,ig)) 

                  if (lay == laysolfr .and. isolvar < 0) &
                     sfluxzen(ig,icol) = sfluxref(ig) 
                  if (lay == laysolfr .and. isolvar >= 0 .and. isolvar <= 2) &
                     ssi(ig,icol) = svar_f * facbrght(ig) + &
                                    svar_s * snsptdrk(ig) + &
                                    svar_i * irradnce(ig)
                  if (lay == laysolfr .and. isolvar == 3) &
                     ssi(ig,icol) = svar_f_bnd(ngb(ig)) * facbrght(ig) + &
                                    svar_s_bnd(ngb(ig)) * snsptdrk(ig) + &
                                    svar_i_bnd(ngb(ig)) * irradnce(ig)
                  taur(lay,ig,icol) = tauray  
               enddo
            end if
         enddo
      enddo

   end subroutine taumol16

   !----------------------------------------------------------------------------
   subroutine taumol17(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng17, ngs16
      use rrsw_kg17, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray, strrat

      layreffr = 30
      strrat = 0.364641 

      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  

      do icol = 1,ncol

         ! Lower atmosphere loop
         do lay = 1,nlayers 
            if (lay <= laytrop(icol)) then
            !do lay = 1,laytrop(icol) 
               speccomb = colh2o(lay,icol) + strrat*colco2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               fac000 = (1. - fs) * fac00(lay,icol) 
               fac010 = (1. - fs) * fac10(lay,icol) 
               fac100 =       fs  * fac00(lay,icol) 
               fac110 =       fs  * fac10(lay,icol) 
               fac001 = (1. - fs) * fac01(lay,icol) 
               fac011 = (1. - fs) * fac11(lay,icol) 
               fac101 =       fs  * fac01(lay,icol) 
               fac111 =       fs  * fac11(lay,icol) 
               ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(17) + js
               ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(17) + js
               inds = indself(lay,icol) 
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng17
                  taug(lay,ngs16+ig,icol) = speccomb * &
                     (fac000 * absa(ind0,   ig) + &
                      fac100 * absa(ind0+1, ig) + &
                      fac010 * absa(ind0+9, ig) + &
                      fac110 * absa(ind0+10,ig) + &
                      fac001 * absa(ind1,   ig) + &
                      fac101 * absa(ind1+1, ig) + &
                      fac011 * absa(ind1+9, ig) + &
                      fac111 * absa(ind1+10,ig)) + &
                     colh2o(lay,icol) * &
                     (selffac(lay,icol) * LIN2_ARG1(selfref,inds,ig,selffrac(lay,icol)) + &
                      forfac (lay,icol) * LIN2_ARG1( forref,indf,ig, forfrac(lay,icol)))
                  taur(lay,ngs16+ig,icol) = tauray
               enddo

            else

               speccomb = colh2o(lay,icol) + strrat*colco2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 4. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               fac000 = (1. - fs) * fac00(lay,icol) 
               fac010 = (1. - fs) * fac10(lay,icol) 
               fac100 =       fs  * fac00(lay,icol) 
               fac110 =       fs  * fac10(lay,icol) 
               fac001 = (1. - fs) * fac01(lay,icol) 
               fac011 = (1. - fs) * fac11(lay,icol) 
               fac101 =       fs  * fac01(lay,icol) 
               fac111 =       fs  * fac11(lay,icol) 
               ind0 = ((jp(lay,icol)-13)*5+(jt (lay,icol)-1))*nspb(17) + js
               ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(17) + js
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng17
                  taug(lay,ngs16+ig,icol) = speccomb * &
                     (fac000 * absb(ind0,  ig) + &
                      fac100 * absb(ind0+1,ig) + &
                      fac010 * absb(ind0+5,ig) + &
                      fac110 * absb(ind0+6,ig) + &
                      fac001 * absb(ind1,  ig) + &
                      fac101 * absb(ind1+1,ig) + &
                      fac011 * absb(ind1+5,ig) + &
                      fac111 * absb(ind1+6,ig)) + &
                     colh2o(lay,icol) * &
                     forfac(lay,icol) * LIN2_ARG1(forref,indf,ig,forfrac(lay,icol)) 
!                 ssa(lay,ngs16+ig) = tauray / taug(lay,ngs16+ig)
                  taur(lay,ngs16+ig,icol) = tauray
               enddo
            endif
         enddo
      enddo
        
      ! Upper atmosphere loop
      do icol = 1,ncol      
         laysolfr = nlayers 
         do lay = 2,nlayers
            if (lay > laytrop(icol)) then 
          
               if ((jp(lay-1,icol) < layreffr) .and. (jp(lay,icol) >= layreffr)) then
                  laysolfr = lay
               end if
          
               if (lay == laysolfr) then
              
                  speccomb = colh2o(lay,icol) + strrat*colco2(lay,icol) 
                  specparm = colh2o(lay,icol) / speccomb 
                  if (specparm >= oneminus) specparm = oneminus
                  specmult = 4. * specparm
                  js = 1 + int(specmult)
                  fs = mod(specmult, 1.)
                  do ig = 1,ng17 
                     if (isolvar < 0) &
                        sfluxzen(ngs16+ig,icol) = LIN2_ARG2(sfluxref,ig,js,fs)
                     if (isolvar >= 0 .and. isolvar <= 2) &
                        ssi(ngs16+ig,icol) = &
                           svar_f * LIN2_ARG2(facbrght,ig,js,fs) + &
                           svar_s * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                           svar_i * LIN2_ARG2(irradnce,ig,js,fs)
                     if (isolvar == 3) &
                        ssi(ngs16+ig,icol) = &
                           svar_f_bnd(ngb(ngs16+ig)) * LIN2_ARG2(facbrght,ig,js,fs) + &
                           svar_s_bnd(ngb(ngs16+ig)) * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                           svar_i_bnd(ngb(ngs16+ig)) * LIN2_ARG2(irradnce,ig,js,fs)
                  end do
               end if
            end if
         enddo
      enddo   

   end subroutine taumol17

   !----------------------------------------------------------------------------
   subroutine taumol18(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng18, ngs17
      use rrsw_kg18, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray, strrat

      strrat = 38.9589
      layreffr = 6

      do icol = 1,ncol
         laysolfr = laytrop(icol)
         do lay = 1,laytrop(icol)
            speccomb = colh2o(lay,icol) + strrat*colch4(lay,icol) 
            specparm = colh2o(lay,icol) / speccomb 
            if (specparm >= oneminus) specparm = oneminus
            specmult = 8. * specparm
            js = 1 + int(specmult)
            fs = mod(specmult, 1.)
            if (jp(lay,icol) < layreffr .and. jp(lay+1,icol) >= layreffr) &
               laysolfr = min(lay+1,laytrop(icol))
            do ig = 1,ng18
               if (lay == laysolfr .and. isolvar < 0) &
                  sfluxzen(ngs17+ig,icol) = &
                     sfluxref(ig,js) + fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
               if (lay == laysolfr .and. isolvar >= 0 .and. isolvar <= 2) &
                  ssi(ngs17+ig,icol) = &
                     svar_f * LIN2_ARG2(facbrght,ig,js,fs) + &
                     svar_s * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                     svar_i * LIN2_ARG2(irradnce,ig,js,fs)
               if (lay == laysolfr .and. isolvar == 3) &
                  ssi(ngs17+ig,icol) = &
                     svar_f_bnd(ngb(ngs17+ig)) * LIN2_ARG2(facbrght,ig,js,fs) + &
                     svar_s_bnd(ngb(ngs17+ig)) * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                     svar_i_bnd(ngb(ngs17+ig)) * LIN2_ARG2(irradnce,ig,js,fs)
            end do
         end do
      end do
      
      do icol = 1,ncol

         do lay = 1,nlayers
            if (lay <= laytrop(icol)) then
            !do lay = 1,laytrop(icol) 
       
               speccomb = colh2o(lay,icol) + strrat*colch4(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               fac000 = (1. - fs) * fac00(lay,icol) 
               fac010 = (1. - fs) * fac10(lay,icol) 
               fac100 =       fs  * fac00(lay,icol) 
               fac110 =       fs  * fac10(lay,icol) 
               fac001 = (1. - fs) * fac01(lay,icol) 
               fac011 = (1. - fs) * fac11(lay,icol) 
               fac101 =       fs  * fac01(lay,icol) 
               fac111 =       fs  * fac11(lay,icol) 
               ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(18) + js
               ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(18) + js
               inds = indself(lay,icol) 
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng18
                  taug(lay,ngs17+ig,icol) = speccomb * &
                     (fac000 * absa(ind0,   ig) + &
                      fac100 * absa(ind0+1, ig) + &
                      fac010 * absa(ind0+9, ig) + &
                      fac110 * absa(ind0+10,ig) + &
                      fac001 * absa(ind1,   ig) + &
                      fac101 * absa(ind1+1, ig) + &
                      fac011 * absa(ind1+9, ig) + &
                      fac111 * absa(ind1+10,ig)) + &
                     colh2o(lay,icol) * &
                     (selffac(lay,icol) * LIN2_ARG1(selfref,inds,ig,selffrac(lay,icol)) + &
                      forfac (lay,icol) * LIN2_ARG1( forref,indf,ig, forfrac(lay,icol)))
!                 ssa(lay,ngs17+ig) = tauray / taug(lay,ngs17+ig)
                  taur(lay,ngs17+ig,icol) = tauray
               enddo

            else

               ! Upper atmosphere loop
               !do lay = laytrop(icol) +1, nlayers
               ind0 = ((jp(lay,icol)-13)*5+(jt (lay,icol)-1))*nspb(18) + 1
               ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(18) + 1
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng18
                  taug(lay,ngs17+ig,icol) = colch4(lay,icol) * &
                     (fac00(lay,icol) * absb(ind0,  ig) + &
                      fac10(lay,icol) * absb(ind0+1,ig) + &
                      fac01(lay,icol) * absb(ind1,  ig) + &	  
                      fac11(lay,icol) * absb(ind1+1,ig)) 
!                 ssa(lay,ngs17+ig) = tauray / taug(lay,ngs17+ig)
                  taur(lay,ngs17+ig,icol) = tauray
               enddo
            end if
         enddo
      enddo
       
   end subroutine taumol18

   !----------------------------------------------------------------------------
   subroutine taumol19(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng19, ngs18
      use rrsw_kg19, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray, strrat

      strrat = 5.49281 
      layreffr = 3      
      
      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  

      do icol = 1,ncol

         laysolfr = laytrop(icol) 
  
         ! Lower atmosphere loop      
         do lay = 1,laytrop(icol) 
            
            if (jp(lay,icol) < layreffr .and. jp(lay+1,icol) >= layreffr) &
               laysolfr = min(lay+1,laytrop(icol))
     
            if (lay == laysolfr) then 
               speccomb = colh2o(lay,icol) + strrat*colco2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
        
               do ig = 1,ng19
                  if (isolvar < 0) &
                     sfluxzen(ngs18+ig,icol) = &
                        sfluxref(ig,js) + fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
                  if (isolvar >= 0 .and. isolvar <= 2) &
                     ssi(ngs18+ig,icol) = &
                        svar_f * LIN2_ARG2(facbrght,ig,js,fs) + &
                        svar_s * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                        svar_i * LIN2_ARG2(irradnce,ig,js,fs)
                  if (isolvar == 3) &
                     ssi(ngs18+ig,icol) = &
                        svar_f_bnd(ngb(ngs18+ig)) * LIN2_ARG2(facbrght,ig,js,fs) + &
                        svar_s_bnd(ngb(ngs18+ig)) * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                        svar_i_bnd(ngb(ngs18+ig)) * LIN2_ARG2(irradnce,ig,js,fs)
               end do
            end if
         end do
      end do
      
      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  

      do icol = 1,ncol

         ! Lower atmosphere loop      
         do lay = 1,nlayers
            if (lay <= laytrop(icol)) then
       
               speccomb = colh2o(lay,icol) + strrat*colco2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               fac000 = (1. - fs) * fac00(lay,icol) 
               fac010 = (1. - fs) * fac10(lay,icol) 
               fac100 =       fs  * fac00(lay,icol) 
               fac110 =       fs  * fac10(lay,icol) 
               fac001 = (1. - fs) * fac01(lay,icol) 
               fac011 = (1. - fs) * fac11(lay,icol) 
               fac101 =       fs  * fac01(lay,icol) 
               fac111 =       fs  * fac11(lay,icol) 
               ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(19) + js
               ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(19) + js
               inds = indself(lay,icol) 
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1 , ng19
                  taug(lay,ngs18+ig,icol) = speccomb * &
                     (fac000 * absa(ind0,   ig) + &
                      fac100 * absa(ind0+1, ig) + &
                      fac010 * absa(ind0+9, ig) + &
                      fac110 * absa(ind0+10,ig) + &
                      fac001 * absa(ind1,   ig) + &
                      fac101 * absa(ind1+1, ig) + &
                      fac011 * absa(ind1+9, ig) + &
                      fac111 * absa(ind1+10,ig)) + &
                     colh2o(lay,icol) * &
                     (selffac(lay,icol) * LIN2_ARG1(selfref,inds,ig,selffrac(lay,icol)) + &
                      forfac (lay,icol) * LIN2_ARG1( forref,indf,ig, forfrac(lay,icol)))
!                 ssa(lay,ngs18+ig) = tauray / taug(lay,ngs18+ig)
                  taur(lay,ngs18+ig,icol) = tauray   
               enddo
            else

               ! Upper atmosphere loop
               ind0 = ((jp(lay,icol)-13)*5+(jt (lay,icol)-1))*nspb(19) + 1
               ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(19) + 1
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng19
                  taug(lay,ngs18+ig,icol) = colco2(lay,icol) * &
                     (fac00(lay,icol) * absb(ind0,  ig) + &
                      fac10(lay,icol) * absb(ind0+1,ig) + &
                      fac01(lay,icol) * absb(ind1,  ig) + &
                      fac11(lay,icol) * absb(ind1+1,ig)) 
!                 ssa(lay,ngs18+ig) = tauray / taug(lay,ngs18+ig) 
                  taur(lay,ngs18+ig,icol) = tauray   
               enddo
            end if
         enddo
      enddo

   end subroutine taumol19

   !----------------------------------------------------------------------------
   subroutine taumol20(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng20, ngs19
      use rrsw_kg20, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, absch4, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      implicit none

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray

      layreffr = 3

      do icol = 1,ncol
         laysolfr = laytrop(icol)
         do lay = 1,laytrop(icol)
            if (jp(lay,icol) < layreffr .and. jp(lay+1,icol) >= layreffr) &
               laysolfr = min(lay+1,laytrop(icol))
            if (lay == laysolfr) then 
               do ig = 1,ng20 
                  if (isolvar < 0) &
                     sfluxzen(ngs19+ig,icol) = sfluxref(ig) 
                  if (isolvar >= 0 .and. isolvar <= 2) &
                     ssi(ngs19+ig,icol) = svar_f * facbrght(ig) + &
                                          svar_s * snsptdrk(ig) + &
                                          svar_i * irradnce(ig)
                  if (isolvar == 3) &
                     ssi(ngs19+ig,icol) = svar_f_bnd(ngb(ngs19+ig)) * facbrght(ig) + &
                                          svar_s_bnd(ngb(ngs19+ig)) * snsptdrk(ig) + &
                                          svar_i_bnd(ngb(ngs19+ig)) * irradnce(ig)
               end do
            end if
         end do
      end do

      do icol = 1,ncol

         do lay = 1,nlayers 
            if (lay <= laytrop(icol)) then
         
               ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(20) + 1
               ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(20) + 1
               inds = indself(lay,icol) 
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng20
                  taug(lay,ngs19+ig,icol) = colh2o(lay,icol) * &
                     ((fac00(lay,icol) * absa(ind0,  ig) + &
                       fac10(lay,icol) * absa(ind0+1,ig) + &
                       fac01(lay,icol) * absa(ind1,  ig) + &
                       fac11(lay,icol) * absa(ind1+1,ig)) + &
                      selffac(lay,icol) * LIN2_ARG1(selfref,inds,ig,selffrac(lay,icol)) + &
                      forfac (lay,icol) * LIN2_ARG1( forref,indf,ig, forfrac(lay,icol))) &
                     + colch4(lay,icol) * absch4(ig)
!                 ssa(lay,ngs19+ig) = tauray / taug(lay,ngs19+ig)
                  taur(lay,ngs19+ig,icol) = tauray 
               enddo
            else

               ! Upper atmosphere loop
               ind0 = ((jp(lay,icol)-13)*5+(jt (lay,icol)-1))*nspb(20) + 1
               ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(20) + 1
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng20
                  taug(lay,ngs19+ig,icol) = colh2o(lay,icol) * &
                     (fac00(lay,icol) * absb(ind0,  ig) + &
                      fac10(lay,icol) * absb(ind0+1,ig) + &
                      fac01(lay,icol) * absb(ind1,  ig) + &
                      fac11(lay,icol) * absb(ind1+1,ig) + &
                      forfac(lay,icol) * LIN2_ARG1(forref,indf,ig,forfrac(lay,icol))) &
                     + colch4(lay,icol) * absch4(ig)
!                 ssa(lay,ngs19+ig) = tauray / taug(lay,ngs19+ig)
                  taur(lay,ngs19+ig,icol) = tauray 
               enddo
            end if
         enddo
      enddo

   end subroutine taumol20

   !----------------------------------------------------------------------------
   subroutine taumol21(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng21, ngs20
      use rrsw_kg21, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray, strrat
        
      strrat = 0.0045321 
      layreffr = 8
        
      do icol = 1,ncol
         laysolfr = laytrop(icol)        
         do lay=1,laytrop(icol)
            if (jp(lay,icol) < layreffr .and. jp(lay+1,icol) >= layreffr) &
               laysolfr = min(lay+1,laytrop(icol))
            if (lay == laysolfr) then 
               speccomb = colh2o(lay,icol) + strrat*colco2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               do ig = 1,ng21
                  if (isolvar < 0) &
                     sfluxzen(ngs20+ig,icol) = &
                        sfluxref(ig,js) + fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
                  if (isolvar >= 0 .and. isolvar <= 2) &
                     ssi(ngs20+ig,icol) = &
                        svar_f * LIN2_ARG2(facbrght,ig,js,fs) + &
                        svar_s * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                        svar_i * LIN2_ARG2(irradnce,ig,js,fs)
                  if (isolvar == 3) &
                     ssi(ngs20+ig,icol) = &
                        svar_f_bnd(ngb(ngs20+ig)) * LIN2_ARG2(facbrght,ig,js,fs) + &
                        svar_s_bnd(ngb(ngs20+ig)) * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                        svar_i_bnd(ngb(ngs20+ig)) * LIN2_ARG2(irradnce,ig,js,fs)
               end do
            end if
         end do
      end do        

      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  

      do icol = 1,ncol  
      
         ! Lower atmosphere loop
         do lay = 1,nlayers 
            if (lay <= laytrop(icol)) then 
        
               speccomb = colh2o(lay,icol) + strrat*colco2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               fac000 = (1. - fs) * fac00(lay,icol) 
               fac010 = (1. - fs) * fac10(lay,icol) 
               fac100 =       fs  * fac00(lay,icol) 
               fac110 =       fs  * fac10(lay,icol) 
               fac001 = (1. - fs) * fac01(lay,icol) 
               fac011 = (1. - fs) * fac11(lay,icol) 
               fac101 =       fs  * fac01(lay,icol) 
               fac111 =       fs  * fac11(lay,icol) 
               ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(21) + js
               ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(21) + js
               inds = indself(lay,icol) 
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng21
                  taug(lay,ngs20+ig,icol) = speccomb * &
                     (fac000 * absa(ind0,   ig) + &
                      fac100 * absa(ind0+1, ig) + &
                      fac010 * absa(ind0+9, ig) + &
                      fac110 * absa(ind0+10,ig) + &
                      fac001 * absa(ind1,   ig) + &
                      fac101 * absa(ind1+1, ig) + &
                      fac011 * absa(ind1+9, ig) + &
                      fac111 * absa(ind1+10,ig)) + &
                     colh2o(lay,icol) * &
                     (selffac(lay,icol) * LIN2_ARG1(selfref,inds,ig,selffrac(lay,icol)) + &
                      forfac (lay,icol) * LIN2_ARG1( forref,indf,ig, forfrac(lay,icol)))
!                 ssa(lay,ngs20+ig) = tauray / taug(lay,ngs20+ig)
                  taur(lay,ngs20+ig,icol) = tauray
               enddo

            else

               ! Upper atmosphere loop
               speccomb = colh2o(lay,icol) + strrat*colco2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 4. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               fac000 = (1. - fs) * fac00(lay,icol) 
               fac010 = (1. - fs) * fac10(lay,icol) 
               fac100 =       fs  * fac00(lay,icol) 
               fac110 =       fs  * fac10(lay,icol) 
               fac001 = (1. - fs) * fac01(lay,icol) 
               fac011 = (1. - fs) * fac11(lay,icol) 
               fac101 =       fs  * fac01(lay,icol) 
               fac111 =       fs  * fac11(lay,icol) 
               ind0 = ((jp(lay,icol) -13)*5+(jt (lay,icol)-1))*nspb(21) + js
               ind1 = ((jp(lay,icol) -12)*5+(jt1(lay,icol)-1))*nspb(21) + js
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng21
                  taug(lay,ngs20+ig,icol) = speccomb * &
                     (fac000 * absb(ind0,  ig) + &
                      fac100 * absb(ind0+1,ig) + &
                      fac010 * absb(ind0+5,ig) + &
                      fac110 * absb(ind0+6,ig) + &
                      fac001 * absb(ind1,  ig) + &
                      fac101 * absb(ind1+1,ig) + &
                      fac011 * absb(ind1+5,ig) + &
                      fac111 * absb(ind1+6,ig)) + &
                     colh2o(lay,icol) * &
                     forfac(lay,icol) * LIN2_ARG1(forref,indf,ig,forfrac(lay,icol)) 
!                 ssa(lay,ngs20+ig) = tauray / taug(lay,ngs20+ig)
                  taur(lay,ngs20+ig,icol) = tauray
               enddo
            end if
         enddo
      enddo

   end subroutine taumol21

   !----------------------------------------------------------------------------
   subroutine taumol22(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng22, ngs21
      use rrsw_kg22, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray, o2adj, o2cont, strrat

      ! The following factor is the ratio of total O2 band intensity (lines 
      ! and Mate continuum) to O2 band intensity (line only).  It is needed
      ! to adjust the optical depths since the k's include only lines.
      o2adj = 1.6
      
      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  

      strrat = 0.022708
      layreffr = 2
      
      do icol = 1,ncol

         laysolfr = laytrop(icol) 

         ! Lower atmosphere loop
         do lay = 1,laytrop(icol) 
            if (jp(lay,icol) < layreffr .and. jp(lay+1,icol) >= layreffr) &
               laysolfr = min(lay+1,laytrop(icol))
                 
            if (lay == laysolfr) then 
               speccomb = colh2o(lay,icol) + o2adj*strrat*colo2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
!              odadj = specparm + o2adj * (1. - specparm)
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               do ig = 1,ng22                                 
                  if (isolvar < 0) &
                     sfluxzen(ngs21+ig,icol) = &
                        sfluxref(ig,js) + fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
                  if (isolvar >= 0 .and. isolvar <= 2) &
                     ssi(ngs21+ig,icol) = &
                        svar_f * LIN2_ARG2(facbrght,ig,js,fs) + &
                        svar_s * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                        svar_i * LIN2_ARG2(irradnce,ig,js,fs)
                  if (isolvar == 3) &
                     ssi(ngs21+ig,icol) = &
                        svar_f_bnd(ngb(ngs21+ig)) * LIN2_ARG2(facbrght,ig,js,fs) + &
                        svar_s_bnd(ngb(ngs21+ig)) * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                        svar_i_bnd(ngb(ngs21+ig)) * LIN2_ARG2(irradnce,ig,js,fs)
               end do
            end if
         end do
      end do
 
      do icol = 1,ncol

         laysolfr = laytrop(icol) 

         ! Lower atmosphere loop

         do lay = 1,nlayers 
            if (lay<=laytrop(icol)) then
  
               o2cont = 4.35e-4 *colo2(lay,icol) /(350.0 *2.0 )
               speccomb = colh2o(lay,icol) + o2adj*strrat*colo2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
!              odadj = specparm + o2adj * (1. - specparm)
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               fac000 = (1. - fs) * fac00(lay,icol) 
               fac010 = (1. - fs) * fac10(lay,icol) 
               fac100 =       fs  * fac00(lay,icol) 
               fac110 =       fs  * fac10(lay,icol) 
               fac001 = (1. - fs) * fac01(lay,icol) 
               fac011 = (1. - fs) * fac11(lay,icol) 
               fac101 =       fs  * fac01(lay,icol) 
               fac111 =       fs  * fac11(lay,icol) 
               ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(22) + js
               ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(22) + js
               inds = indself(lay,icol) 
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng22
                  taug(lay,ngs21+ig,icol) = speccomb * &
                     (fac000 * absa(ind0,   ig) + &
                      fac100 * absa(ind0+1, ig) + &
                      fac010 * absa(ind0+9, ig) + &
                      fac110 * absa(ind0+10,ig) + &
                      fac001 * absa(ind1,   ig) + &
                      fac101 * absa(ind1+1, ig) + &
                      fac011 * absa(ind1+9, ig) + &
                      fac111 * absa(ind1+10,ig)) + &
                     colh2o(lay,icol) * &
                     (selffac(lay,icol) * LIN2_ARG1(selfref,inds,ig,selffrac(lay,icol)) + &
                      forfac (lay,icol) * LIN2_ARG1( forref,indf,ig, forfrac(lay,icol))) &
                     + o2cont
!                 ssa(lay,ngs21+ig) = tauray / taug(lay,ngs21+ig)
                  taur(lay,ngs21+ig,icol) = tauray
               enddo

            else

               ! Upper atmosphere loop
               o2cont = 4.35e-4 *colo2(lay,icol) /(350.0 *2.0 )
               ind0 = ((jp(lay,icol)-13)*5+(jt (lay,icol)-1))*nspb(22) + 1
               ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(22) + 1
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng22
                  taug(lay,ngs21+ig,icol) = colo2(lay,icol) * o2adj * &
                     (fac00(lay,icol) * absb(ind0,  ig) + &
                      fac10(lay,icol) * absb(ind0+1,ig) + &
                      fac01(lay,icol) * absb(ind1,  ig) + &
                      fac11(lay,icol) * absb(ind1+1,ig)) + &
                     o2cont
!                 ssa(lay,ngs21+ig) = tauray / taug(lay,ngs21+ig)
                  taur(lay,ngs21+ig,icol) = tauray
               enddo
            end if
         enddo
      enddo

   end subroutine taumol22

   !----------------------------------------------------------------------------
   subroutine taumol23(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 23:  8050-12850 cm-1 (low - h2o; high - nothing)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng23, ngs22
      use rrsw_kg23, only : absa, ka, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray, givfac

      ! Average Giver et al. correction factor for this band.
      givfac = 1.029

      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  

      layreffr = 6
      
      do icol = 1,ncol

         laysolfr = laytrop(icol) 

         ! Lower atmosphere loop
         do lay = 1,laytrop(icol) 
            if (jp(lay,icol) < layreffr .and. jp(lay+1,icol) >= layreffr) &
               laysolfr = min(lay+1,laytrop(icol))

            if (lay == laysolfr) then 
               do ig = 1,ng23
                  if (isolvar < 0) &
                     sfluxzen(ngs22+ig,icol) = sfluxref(ig) 
                  if (isolvar >= 0 .and. isolvar <= 2) &
                     ssi(ngs22+ig,icol) = svar_f * facbrght(ig) + &
                                          svar_s * snsptdrk(ig) + &
                                          svar_i * irradnce(ig)
                  if (isolvar == 3) &
                     ssi(ngs22+ig,icol) = svar_f_bnd(ngb(ngs22+ig)) * facbrght(ig) + &
                                          svar_s_bnd(ngb(ngs22+ig)) * snsptdrk(ig) + &
                                          svar_i_bnd(ngb(ngs22+ig)) * irradnce(ig)
               end do
            end if
         end do
      end do      
      
      do icol = 1,ncol

         ! Lower atmosphere loop
         do lay = 1,nlayers 
            if (lay <= laytrop(icol)) then
               if (jp(lay,icol) < layreffr .and. jp(lay+1,icol) >= layreffr) &
                  laysolfr = min(lay+1,laytrop(icol))
               ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(23) + 1
               ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(23) + 1
               inds = indself(lay,icol) 
               indf = indfor(lay,icol) 

               do ig = 1,ng23
                  tauray = colmol(lay,icol) * rayl(ig)
                  taug(lay,ngs22+ig,icol) = colh2o(lay,icol) * (givfac * &
                     (fac00(lay,icol) * absa(ind0,  ig) + &
                      fac10(lay,icol) * absa(ind0+1,ig) + &
                      fac01(lay,icol) * absa(ind1,  ig) + &
                      fac11(lay,icol) * absa(ind1+1,ig)) + &
                     selffac(lay,icol) * LIN2_ARG1(selfref,inds,ig,selffrac(lay,icol)) + &
                     forfac (lay,icol) * LIN2_ARG1( forref,indf,ig, forfrac(lay,icol)))
!                 ssa(lay,ngs22+ig) = tauray / taug(lay,ngs22+ig)
                  taur(lay,ngs22+ig,icol) = tauray
               enddo
            else

               ! Upper atmosphere loop
               do ig = 1,ng23
!                 taug(lay,ngs22+ig) = colmol(lay) * rayl(ig)
!                 ssa(lay,ngs22+ig) = 1.
                  taug(lay,ngs22+ig,icol) = 0. 
                  taur(lay,ngs22+ig,icol) = colmol(lay,icol) * rayl(ig) 
               enddo
            end if
         enddo
      enddo

   end subroutine taumol23

   !----------------------------------------------------------------------------
   subroutine taumol24(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng24, ngs23
      use rrsw_kg24, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, abso3a, abso3b, rayla, raylb, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray, strrat

      strrat = 0.124692 
      layreffr = 1   
        
      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  

      do icol = 1,ncol

         laysolfr = laytrop(icol) 

         ! Lower atmosphere loop
         do lay = 1,laytrop(icol) 
            if (jp(lay,icol) < layreffr .and. jp(lay+1,icol) >= layreffr) &
               laysolfr = min(lay+1,laytrop(icol))
            if (lay == laysolfr) then
               speccomb = colh2o(lay,icol) + strrat*colo2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               do ig = 1,ng24
                  if (isolvar < 0) &
                  sfluxzen(ngs23+ig,icol) = sfluxref(ig,js) + &
                     fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
                  if (isolvar >= 0 .and. isolvar <= 2) &
                     ssi(ngs23+ig,icol) = &
                        svar_f * LIN2_ARG2(facbrght,ig,js,fs) + &
                        svar_s * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                        svar_i * LIN2_ARG2(irradnce,ig,js,fs)
                  if (isolvar == 3) &
                     ssi(ngs23+ig,icol) = &
                        svar_f_bnd(ngb(ngs23+ig)) * LIN2_ARG2(facbrght,ig,js,fs) + &
                        svar_s_bnd(ngb(ngs23+ig)) * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                        svar_i_bnd(ngb(ngs23+ig)) * LIN2_ARG2(irradnce,ig,js,fs)
               end do
            end if
         end do
      end do
     
      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  

      do icol = 1,ncol

         ! Lower atmosphere loop
         do lay = 1,nlayers 
            if (lay <= laytrop(icol)) then

               speccomb = colh2o(lay,icol) + strrat*colo2(lay,icol) 
               specparm = colh2o(lay,icol) / speccomb 
               if (specparm >= oneminus) specparm = oneminus
               specmult = 8. * specparm
               js = 1 + int(specmult)
               fs = mod(specmult, 1.)
               fac000 = (1. - fs) * fac00(lay,icol) 
               fac010 = (1. - fs) * fac10(lay,icol) 
               fac100 =       fs  * fac00(lay,icol) 
               fac110 =       fs  * fac10(lay,icol) 
               fac001 = (1. - fs) * fac01(lay,icol) 
               fac011 = (1. - fs) * fac11(lay,icol) 
               fac101 =       fs  * fac01(lay,icol) 
               fac111 =       fs  * fac11(lay,icol) 
               ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(24) + js
               ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(24) + js
               inds = indself(lay,icol) 
               indf = indfor (lay,icol) 

               do ig = 1,ng24
                  tauray = colmol(lay,icol) * (rayla(ig,js) + fs * (rayla(ig,js+1) - rayla(ig,js)))
                  taug(lay,ngs23+ig,icol) = speccomb * &
                     (fac000 * absa(ind0,   ig) + &
                      fac100 * absa(ind0+1, ig) + &
                      fac010 * absa(ind0+9, ig) + &
                      fac110 * absa(ind0+10,ig) + &
                      fac001 * absa(ind1,   ig) + &
                      fac101 * absa(ind1+1, ig) + &
                      fac011 * absa(ind1+9, ig) + &
                      fac111 * absa(ind1+10,ig)) + &
                     colo3(lay,icol) * abso3a(ig) + &
                     colh2o(lay,icol) * & 
                     (selffac(lay,icol) * LIN2_ARG1(selfref,inds,ig,selffrac(lay,icol)) + &
                      forfac (lay,icol) * LIN2_ARG1( forref,indf,ig, forfrac(lay,icol)))
!                 ssa(lay,ngs23+ig) = tauray / taug(lay,ngs23+ig)
                  taur(lay,ngs23+ig,icol) = tauray
               enddo

            else

               ! Upper atmosphere loop
               ind0 = ((jp(lay,icol)-13)*5+(jt (lay,icol)-1))*nspb(24) + 1
               ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(24) + 1

               do ig = 1,ng24
                  tauray = colmol(lay,icol) * raylb(ig)
                  taug(lay,ngs23+ig,icol) = colo2(lay,icol) * &
                     (fac00(lay,icol) * absb(ind0,  ig) + &
                      fac10(lay,icol) * absb(ind0+1,ig) + &
                      fac01(lay,icol) * absb(ind1,  ig) + &
                      fac11(lay,icol) * absb(ind1+1,ig)) + &
                     colo3(lay,icol) * abso3b(ig)
!                 ssa(lay,ngs23+ig) = tauray / taug(lay,ngs23+ig)
                  taur(lay,ngs23+ig,icol) = tauray
               enddo
            endif
         enddo
      enddo

   end subroutine taumol24

   !----------------------------------------------------------------------------
   subroutine taumol25(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 25:  16000-22650 cm-1 (low - h2o; high - nothing)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng25, ngs24
      use rrsw_kg25, only : absa, ka, &
                            sfluxref, abso3a, abso3b, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray

      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  
       
      do icol = 1,ncol

         layreffr = 2
         laysolfr = laytrop(icol) 

         ! Lower atmosphere loop
         do lay = 1,laytrop(icol) 
            if (jp(lay,icol) < layreffr .and. jp(lay+1,icol) >= layreffr) &
               laysolfr = min(lay+1,laytrop(icol))
            ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(25) + 1
            ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(25) + 1

            do ig = 1,ng25
               tauray = colmol(lay,icol) * rayl(ig)
               taug(lay,ngs24+ig,icol) = colh2o(lay,icol) * &
                  (fac00(lay,icol) * absa(ind0,  ig) + &
                   fac10(lay,icol) * absa(ind0+1,ig) + &
                   fac01(lay,icol) * absa(ind1,  ig) + &
                   fac11(lay,icol) * absa(ind1+1,ig)) + &
                  colo3(lay,icol) * abso3a(ig) 
!              ssa(lay,ngs24+ig) = tauray / taug(lay,ngs24+ig)
               if (lay == laysolfr .and. isolvar < 0) &
                  sfluxzen(ngs24+ig,icol) = sfluxref(ig) 
               if (lay == laysolfr .and. isolvar >= 0 .and. isolvar <= 2) &
                  ssi(ngs24+ig,icol) = svar_f * facbrght(ig) + &
                                       svar_s * snsptdrk(ig) + &
                                       svar_i * irradnce(ig)
               if (lay == laysolfr .and. isolvar == 3) &
                  ssi(ngs24+ig,icol) = svar_f_bnd(ngb(ngs24+ig)) * facbrght(ig) + &
                                       svar_s_bnd(ngb(ngs24+ig)) * snsptdrk(ig) + &
                                       svar_i_bnd(ngb(ngs24+ig)) * irradnce(ig)
               taur(lay,ngs24+ig,icol) = tauray
            enddo
         enddo

         ! Upper atmosphere loop
         do lay = laytrop(icol) +1, nlayers
            do ig = 1,ng25
               tauray = colmol(lay,icol) * rayl(ig)
               taug(lay,ngs24+ig,icol) = colo3(lay,icol) * abso3b(ig) 
!              ssa(lay,ngs24+ig) = tauray / taug(lay,ngs24+ig)
               taur(lay,ngs24+ig,icol) = tauray
            enddo
         enddo
      enddo

   end subroutine taumol25

   !----------------------------------------------------------------------------
   subroutine taumol26(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng26, ngs25
      use rrsw_kg26, only : sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray

      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  
       
      do icol = 1,ncol

         laysolfr = laytrop(icol) 

         ! Lower atmosphere loop
         do lay = 1,laytrop(icol) 
            do ig = 1,ng26 
!              taug(lay,ngs25+ig) = colmol(lay) * rayl(ig)
!              ssa(lay,ngs25+ig) = 1.
               if (lay == laysolfr .and. isolvar < 0) &
                  sfluxzen(ngs25+ig,icol) = sfluxref(ig) 
               if (lay == laysolfr .and. isolvar >= 0 .and. isolvar <= 2) &
                  ssi(ngs25+ig,icol) = svar_f * facbrght(ig) + &
                                       svar_s * snsptdrk(ig) + &
                                       svar_i * irradnce(ig)
               if (lay == laysolfr .and. isolvar == 3) &
                  ssi(ngs25+ig,icol) = svar_f_bnd(ngb(ngs25+ig)) * facbrght(ig) + &
                                       svar_s_bnd(ngb(ngs25+ig)) * snsptdrk(ig) + &
                                       svar_i_bnd(ngb(ngs25+ig)) * irradnce(ig)
               taug(lay,ngs25+ig,icol) = 0. 
               taur(lay,ngs25+ig,icol) = colmol(lay,icol) * rayl(ig) 
            enddo
         enddo

         ! Upper atmosphere loop
         do lay = laytrop(icol)+1, nlayers
            do ig = 1,ng26
!              taug(lay,ngs25+ig) = colmol(lay) * rayl(ig)
!              ssa(lay,ngs25+ig) = 1.
               taug(lay,ngs25+ig,icol) = 0. 
               taur(lay,ngs25+ig,icol) = colmol(lay,icol) * rayl(ig) 
            enddo
         enddo
      enddo

   end subroutine taumol26

   !----------------------------------------------------------------------------
   subroutine taumol27(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 27:  29000-38000 cm-1 (low - o3; high - o3)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng27, ngs26
      use rrsw_kg27, only : absa, ka, absb, kb, sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray, scalekur

      do icol = 1,ncol

         ! Kurucz solar source function
         ! The values in sfluxref were obtained using the "low resolution"
         ! version of the Kurucz solar source function.  For unknown reasons,
         ! the total irradiance in this band differs from the corresponding
         ! total in the "high-resolution" version of the Kurucz function.
         ! Therefore, these values are scaled below by the factor SCALEKUR.

         scalekur = 50.15/48.37

         ! Compute the optical depth by interpolating in ln(pressure), 
         ! temperature, and appropriate species.  Below LAYTROP, the water
         ! vapor self-continuum is interpolated (in temperature) separately.  

         layreffr = 32

         ! Lower atmosphere loop
         do lay = 1,laytrop(icol) 
            ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(27) + 1
            ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(27) + 1

            do ig = 1,ng27
               tauray = colmol(lay,icol) * rayl(ig)
               taug(lay,ngs26+ig,icol) = colo3(lay,icol) * &
                  (fac00(lay,icol) * absa(ind0,  ig) + &
                   fac10(lay,icol) * absa(ind0+1,ig) + &
                   fac01(lay,icol) * absa(ind1,  ig) + &
                   fac11(lay,icol) * absa(ind1+1,ig))
!              ssa(lay,ngs26+ig) = tauray / taug(lay,ngs26+ig)
               taur(lay,ngs26+ig,icol) = tauray
            enddo
         enddo

         laysolfr = nlayers

         ! Upper atmosphere loop
         do lay = laytrop(icol)+1, nlayers
            if (jp(lay-1,icol) < layreffr .and. jp(lay,icol) >= layreffr) &
               laysolfr = lay
            ind0 = ((jp(lay,icol)-13)*5+(jt (lay,icol)-1))*nspb(27) + 1
            ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(27) + 1

            do ig = 1,ng27
               tauray = colmol(lay,icol) * rayl(ig)
               taug(lay,ngs26+ig,icol) = colo3(lay,icol) * &
                  (fac00(lay,icol) * absb(ind0,  ig) + &
                   fac10(lay,icol) * absb(ind0+1,ig) + &
                   fac01(lay,icol) * absb(ind1,  ig) + & 
                   fac11(lay,icol) * absb(ind1+1,ig))
!              ssa(lay,ngs26+ig) = tauray / taug(lay,ngs26+ig)
               if (lay == laysolfr .and. isolvar < 0) &
                  sfluxzen(ngs26+ig,icol) = scalekur * sfluxref(ig) 
               if (lay == laysolfr .and. isolvar >= 0 .and. isolvar <= 2) &
                  ssi(ngs26+ig,icol) = svar_f * facbrght(ig) + &
                                       svar_s * snsptdrk(ig) + &
                                       svar_i * irradnce(ig)
               if (lay == laysolfr .and. isolvar == 3) &
                  ssi(ngs26+ig,icol) = svar_f_bnd(ngb(ngs26+ig)) * facbrght(ig) + &
                                       svar_s_bnd(ngb(ngs26+ig)) * snsptdrk(ig) + &
                                       svar_i_bnd(ngb(ngs26+ig)) * irradnce(ig)
               taur(lay,ngs26+ig,icol) = tauray
            enddo
         enddo
      enddo

   end subroutine taumol27

   !----------------------------------------------------------------------------
   subroutine taumol28(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 28:  38000-50000 cm-1 (low - o3,o2; high - o3,o2)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng28, ngs27
      use rrsw_kg28, only : absa, ka, absb, kb, sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray, strrat

      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  
       
      do icol = 1,ncol

         strrat = 6.67029e-07

         ! Per Eli:
         !  Robin Hogan recently sent me some results that showed that the
         !  stratopause heating rates in RRTMG in the 38000-5000 cm-1 band are
         !  negatively biased wrt to LBL calculations in that same region, so
         !  invoking an additional spectral region for the issue you are seeing
         !  may not be necessary.  He hacked the code to get the correct answer in
         !  this band by removing the FS interpolation in taumol for band 28 and
         !  replacing it by always using the values in the 5th (i.e. last) location
         !  in the solar irradiance array. My analysis suggested a different (and
         !  better, assuming it works) remedy, which he hasn’t yet commented on or
         !  tried. (note:  I sent him this suggestion just yesterday.) I think I may
         !  have set the reference solar mapping layer (LAYREFFR) too high in that
         !  band. Instead of 58, I think I should have set it around 40. I suggest
         !  you try that and see if it gets rid of the bias.

         ! Followup from Eli:
         !  RRTMG_SW v4.10 has this set to 42

         layreffr = 42

         ! Lower atmosphere loop
         do lay = 1,laytrop(icol) 
            speccomb = colo3(lay,icol) + strrat*colo2(lay,icol) 
            specparm = colo3(lay,icol) / speccomb 
            if (specparm >= oneminus) specparm = oneminus
            specmult = 8. * specparm
            js = 1 + int(specmult)
            fs = mod(specmult, 1.)
            fac000 = (1. - fs) * fac00(lay,icol) 
            fac010 = (1. - fs) * fac10(lay,icol) 
            fac100 =       fs  * fac00(lay,icol) 
            fac110 =       fs  * fac10(lay,icol) 
            fac001 = (1. - fs) * fac01(lay,icol) 
            fac011 = (1. - fs) * fac11(lay,icol) 
            fac101 =       fs  * fac01(lay,icol) 
            fac111 =       fs  * fac11(lay,icol) 
            ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(28) + js
            ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(28) + js
            tauray = colmol(lay,icol) * rayl

            do ig = 1,ng28
               taug(lay,ngs27+ig,icol) = speccomb * &
                  (fac000 * absa(ind0,   ig) + &
                   fac100 * absa(ind0+1, ig) + &
                   fac010 * absa(ind0+9, ig) + &
                   fac110 * absa(ind0+10,ig) + &
                   fac001 * absa(ind1,   ig) + &
                   fac101 * absa(ind1+1, ig) + &
                   fac011 * absa(ind1+9, ig) + &
                   fac111 * absa(ind1+10,ig)) 
!              ssa(lay,ngs27+ig) = tauray / taug(lay,ngs27+ig)
               taur(lay,ngs27+ig,icol) = tauray
            enddo
         enddo

         laysolfr = nlayers

         ! Upper atmosphere loop
         do lay = laytrop(icol) +1, nlayers
            if (jp(lay-1,icol) < layreffr .and. jp(lay,icol) >= layreffr) &
               laysolfr = lay
            speccomb = colo3(lay,icol) + strrat*colo2(lay,icol) 
            specparm = colo3(lay,icol) / speccomb 
            if (specparm >= oneminus) specparm = oneminus
            specmult = 4. * specparm
            js = 1 + int(specmult)
            fs = mod(specmult, 1.)
            fac000 = (1. - fs) * fac00(lay,icol) 
            fac010 = (1. - fs) * fac10(lay,icol) 
            fac100 =       fs  * fac00(lay,icol) 
            fac110 =       fs  * fac10(lay,icol) 
            fac001 = (1. - fs) * fac01(lay,icol) 
            fac011 = (1. - fs) * fac11(lay,icol) 
            fac101 =       fs  * fac01(lay,icol) 
            fac111 =       fs  * fac11(lay,icol) 
            ind0 = ((jp(lay,icol)-13)*5+(jt (lay,icol)-1))*nspb(28) + js
            ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(28) + js
            tauray = colmol(lay,icol) * rayl

            do ig = 1,ng28
               taug(lay,ngs27+ig,icol) = speccomb * &
                  (fac000 * absb(ind0,  ig) + &
                   fac100 * absb(ind0+1,ig) + &
                   fac010 * absb(ind0+5,ig) + &
                   fac110 * absb(ind0+6,ig) + &
                   fac001 * absb(ind1,  ig) + &
                   fac101 * absb(ind1+1,ig) + &
                   fac011 * absb(ind1+5,ig) + &
                   fac111 * absb(ind1+6,ig)) 
!              ssa(lay,ngs27+ig) = tauray / taug(lay,ngs27+ig)
               if (lay == laysolfr .and. isolvar < 0) &
                  sfluxzen(ngs27+ig,icol) = &
                     sfluxref(ig,js) + fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
               if (lay == laysolfr .and. isolvar >= 0 .and. isolvar <= 2) &
                  ssi(ngs27+ig,icol) = &
                     svar_f * LIN2_ARG2(facbrght,ig,js,fs) + &
                     svar_s * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                     svar_i * LIN2_ARG2(irradnce,ig,js,fs)
               if (lay == laysolfr .and. isolvar == 3) &
                  ssi(ngs27+ig,icol) = &
                     svar_f_bnd(ngb(ngs27+ig)) * LIN2_ARG2(facbrght,ig,js,fs) + &
                     svar_s_bnd(ngb(ngs27+ig)) * LIN2_ARG2(snsptdrk,ig,js,fs) + &
                     svar_i_bnd(ngb(ngs27+ig)) * LIN2_ARG2(irradnce,ig,js,fs)
               taur(lay,ngs27+ig,icol) = tauray
            enddo
         enddo
      enddo

   end subroutine taumol28

   !----------------------------------------------------------------------------
   subroutine taumol29(pncol, ncol, nlayers, &
                       colh2o, colco2, colch4, colo2, colo3, colmol, &
                       laytrop, jp, jt, jt1, &
                       fac00, fac01, fac10, fac11, &
                       isolvar, svar_f, svar_s, svar_i, &
                       svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                       selffac, selffrac, indself, forfac, forfrac, indfor, &
                       ssi, sfluxzen, taug, taur)
   !----------------------------------------------------------------------------
   !
   !     band 29:  820-2600 cm-1 (low - h2o; high - co2)
   !
   !----------------------------------------------------------------------------

      use parrrsw, only : ng29, ngs28
      use rrsw_kg29, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, absh2o, absco2, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

      integer, intent(in)  :: pncol           ! Dimensioned num of gridcols
      integer, intent(in)  :: ncol            ! Actual number of gridcols
      integer, intent(in)  :: nlayers         ! total number of layers

      integer, intent(in)  :: laytrop     (pncol)            ! tropopause layer index
      integer, intent(in)  :: jp  (nlayers,pncol)
      integer, intent(in)  :: jt  (nlayers,pncol)
      integer, intent(in)  :: jt1 (nlayers,pncol)

      real,    intent(in)  :: colh2o (nlayers,pncol)    ! column amount (h2o)
      real,    intent(in)  :: colco2 (nlayers,pncol)    ! column amount (co2)
      real,    intent(in)  :: colo3  (nlayers,pncol)    ! column amount (o3)
      real,    intent(in)  :: colch4 (nlayers,pncol)    ! column amount (ch4)
      real,    intent(in)  :: colo2  (nlayers,pncol)    ! column amount (o2)
      real,    intent(in)  :: colmol (nlayers,pncol)    ! 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,pncol)
      integer, intent(in) :: indfor   (nlayers,pncol)
      real,    intent(in) :: selffac  (nlayers,pncol)
      real,    intent(in) :: selffrac (nlayers,pncol)
      real,    intent(in) :: forfac   (nlayers,pncol)
      real,    intent(in) :: forfrac  (nlayers,pncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,pncol) &
         :: fac00, fac01, fac10, fac11

      ! Solar variability
      integer, intent(in)  :: isolvar               ! Flag for solar var method
      real,    intent(in)  :: svar_f                ! facular  multiplier
      real,    intent(in)  :: svar_s                ! sunspot  multiplier
      real,    intent(in)  :: svar_i                ! baseline multiplier
      real,    intent(in)  :: svar_f_bnd(jpband)    ! facular  multiplier (by band)
      real,    intent(in)  :: svar_s_bnd(jpband)    ! sunspot  multiplier (by band)
      real,    intent(in)  :: svar_i_bnd(jpband)    ! baseline multiplier (by band)

      real,    intent(out) :: ssi(ngptsw,pncol)           ! spectral solar intensity with solar var
      real,    intent(out) :: sfluxzen(ngptsw,pncol)      ! solar source function
      real,    intent(out) :: taug(nlayers,ngptsw,pncol)  ! Gaseous optical depth 
      real,    intent(out) :: taur(nlayers,ngptsw,pncol)  ! Rayleigh 

      ! ----- Locals -----

      integer :: icol, ig, ind0, ind1, inds, indf, js, lay, laysolfr, layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, fac110, fac111, &
              fs, speccomb, specmult, specparm, tauray

      layreffr = 49  
        
      do icol = 1,ncol
        
         laysolfr = nlayers
         do lay = laytrop(icol) +1, nlayers
            if (jp(lay-1,icol) < layreffr .and. jp(lay,icol) >= layreffr) &
               laysolfr = lay

            if (lay == laysolfr) then 
               do ig = 1,ng29
                  if (isolvar < 0) &
                     sfluxzen(ngs28+ig,icol) = sfluxref(ig) 
                  if (isolvar >= 0 .and. isolvar <= 2) &
                     ssi(ngs28+ig,icol) = svar_f * facbrght(ig) + &
                                          svar_s * snsptdrk(ig) + &
                                          svar_i * irradnce(ig)
                  if (isolvar == 3) &
                     ssi(ngs28+ig,icol) = svar_f_bnd(ngb(ngs28+ig)) * facbrght(ig) + &
                                          svar_s_bnd(ngb(ngs28+ig)) * snsptdrk(ig) + &
                                          svar_i_bnd(ngb(ngs28+ig)) * irradnce(ig)
                end do
            end if
         end do
      end do

      ! Compute the optical depth by interpolating in ln(pressure), 
      ! temperature, and appropriate species.  Below LAYTROP, the water
      ! vapor self-continuum is interpolated (in temperature) separately.  

      do icol = 1,ncol

         ! Lower atmosphere loop
         do lay = 1,nlayers 
            if (lay <= laytrop(icol)) then
               ind0 = ((jp(lay,icol)-1)*5+(jt (lay,icol)-1))*nspa(29) + 1
               ind1 = ( jp(lay,icol)   *5+(jt1(lay,icol)-1))*nspa(29) + 1
               inds = indself(lay,icol) 
               indf = indfor(lay,icol) 
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng29
                  taug(lay,ngs28+ig,icol) = colh2o(lay,icol) * &
                     ((fac00(lay,icol) * absa(ind0,  ig) + &
                       fac10(lay,icol) * absa(ind0+1,ig) + &
                       fac01(lay,icol) * absa(ind1,  ig) + &
                       fac11(lay,icol) * absa(ind1+1,ig)) + &
                      selffac(lay,icol) * LIN2_ARG1(selfref,inds,ig,selffrac(lay,icol)) + &
                      forfac (lay,icol) * LIN2_ARG1( forref,indf,ig, forfrac(lay,icol))) &
                     + colco2(lay,icol) * absco2(ig) 
!                 ssa(lay,ngs28+ig) = tauray / taug(lay,ngs28+ig)
                  taur(lay,ngs28+ig,icol) = tauray
               enddo

            else 

               ! Upper atmosphere loop
               ind0 = ((jp(lay,icol)-13)*5+(jt (lay,icol)-1))*nspb(29) + 1
               ind1 = ((jp(lay,icol)-12)*5+(jt1(lay,icol)-1))*nspb(29) + 1
               tauray = colmol(lay,icol) * rayl

               do ig = 1,ng29
                  taug(lay,ngs28+ig,icol) = colco2(lay,icol) * &
                     (fac00(lay,icol) * absb(ind0,  ig) + &
                      fac10(lay,icol) * absb(ind0+1,ig) + &
                      fac01(lay,icol) * absb(ind1,  ig) + &
                      fac11(lay,icol) * absb(ind1+1,ig)) &  
                     + colh2o(lay,icol) * absh2o(ig) 
!                 ssa(lay,ngs28+ig) = tauray / taug(lay,ngs28+ig)
                  taur(lay,ngs28+ig,icol) = tauray
               enddo
            end if
         enddo
      enddo

   end subroutine taumol29

end module rrtmg_sw_taumol
