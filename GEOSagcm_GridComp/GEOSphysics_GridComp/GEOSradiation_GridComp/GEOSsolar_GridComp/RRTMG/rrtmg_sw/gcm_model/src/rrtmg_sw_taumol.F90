!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

#ifdef _CUDA
#define gpu_device ,device
#else
#define gpu_device 
#endif

      module rrtmg_sw_taumol

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


      use rrsw_con, only: oneminus
      use rrsw_wvn, only: nspa, nspb
      use rrsw_vsn, only: hvrtau, hnamtau

      implicit none

      contains

!----------------------------------------------------------------------------
      subroutine taumol_sw(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           ssi, sfluxzen, taug, taur)

      integer , intent(in) :: ncol
      integer, intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real, intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real , intent(inout) gpu_device :: ssi(:,:)           ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real , intent(inout) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real , intent(inout) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real , intent(inout) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)

! Calculate gaseous optical depth and planck fractions for each spectral band.

      call taumol16(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol17(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol18(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol19(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol20(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol21(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol22(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol23(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol24(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol25(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol26(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol27(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol28(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)

      call taumol29(ncol, nlayers, &
                           colh2o, colco2, colch4, colo2, colo3, colmol, &
                           laytrop, jp, jt, jt1, &
                           fac00, fac01, fac10, fac11, &
                           isolvar, svar_f, svar_s, svar_i, &
                           svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                           selffac, selffrac, indself, forfac, forfrac, indfor, &
                           ssi, sfluxzen, taug, taur)



end subroutine


!----------------------------------------------------------------------------
      subroutine taumol16(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng16
      use rrsw_kg16, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray, strrat1
       integer :: iplon
      strrat1 = 252.131 
      layreffr = 18

!$acc kernels 
 do iplon=1,ncol
! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  



! Lower atmosphere loop
      do lay = 1, nlayers
         if (lay <= laytrop(iplon)) then 
         speccomb = colh2o(iplon,lay)  + strrat1*colch4(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8.*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(16) + js
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(16) + js
         inds = indself(iplon,lay) 
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng16
            taug(iplon,lay,ig)  = speccomb * &
                (fac000 * absa(ind0   ,ig) + &
                 fac100 * absa(ind0 +1,ig) + &
                 fac010 * absa(ind0 +9,ig) + &
                 fac110 * absa(ind0+10,ig) + &
                 fac001 * absa(ind1   ,ig) + &
                 fac101 * absa(ind1 +1,ig) + &
                 fac011 * absa(ind1 +9,ig) + &
                 fac111 * absa(ind1+10,ig)) + &
                 colh2o(iplon,lay)  * &
                 (selffac(iplon,lay)  * (selfref(inds,ig) + &
                 selffrac(iplon,lay)  * &
                 (selfref(inds+1,ig) - selfref(inds,ig))) + &
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig)))) 
!            ssa(lay,ig) = tauray/taug(lay,ig)
            taur(iplon,lay,ig)  = tauray
    
         enddo
         end if
      enddo

     
 end do
 !$acc end kernels

! Upper atmosphere loop
!$acc kernels 
 do iplon=1,ncol
       laysolfr = nlayers
      do lay = 1, nlayers
      if (lay > laytrop(iplon)) then
          !do lay = laytrop(iplon) +1, nlayers
         if (jp(iplon,lay-1)  .lt. layreffr .and. jp(iplon,lay)  .ge. layreffr) then
            laysolfr = lay
         end if
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(16) + 1
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(16) + 1
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng16
            taug(iplon,lay,ig)  = colch4(iplon,lay)  * &
                (fac00(iplon,lay)  * absb(ind0  ,ig) + &
                 fac10(iplon,lay)  * absb(ind0+1,ig) + &
                 fac01(iplon,lay)  * absb(ind1  ,ig) + &
                 fac11(iplon,lay)  * absb(ind1+1,ig)) 

            if (lay .eq. laysolfr .and. isolvar .lt. 0) &
               sfluxzen(iplon,ig) = sfluxref(ig) 
            if (lay .eq. laysolfr .and. isolvar .ge. 0 .and. isolvar .le. 2) &
               ssi(iplon,ig) = svar_f * facbrght(ig) + &
                         svar_s * snsptdrk(ig) + &
                         svar_i * irradnce(ig)
            if (lay .eq. laysolfr .and. isolvar .eq. 3) &
               ssi(iplon,ig) = svar_f_bnd(ngb(ig)) * facbrght(ig) + &
                         svar_s_bnd(ngb(ig)) * snsptdrk(ig) + &
                         svar_i_bnd(ngb(ig)) * irradnce(ig)
            taur(iplon,lay,ig)  = tauray  
         enddo
      end if
         enddo
      enddo

!$acc end kernels
      end subroutine taumol16

!----------------------------------------------------------------------------
      subroutine taumol17(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng17, ngs16
      use rrsw_kg17, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray, strrat
    integer :: iplon

     

   
      layreffr = 30
       strrat = 0.364641 

!$acc kernels loop
 do iplon=1,ncol
! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  

! Lower atmosphere loop
!$acc loop private(js, fs)
      do lay = 1, nlayers 
      if (lay <= laytrop(iplon)) then
          !do lay = 1, laytrop(iplon) 
         speccomb = colh2o(iplon,lay)  + strrat*colco2(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8.*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(17) + js
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(17) + js
         inds = indself(iplon,lay) 
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng17
            taug(iplon,lay,ngs16+ig)  = speccomb * &
                (fac000 * absa(ind0,ig) + &
                 fac100 * absa(ind0+1,ig) + &
                 fac010 * absa(ind0+9,ig) + &
                 fac110 * absa(ind0+10,ig) + &
                 fac001 * absa(ind1,ig) + &
                 fac101 * absa(ind1+1,ig) + &
                 fac011 * absa(ind1+9,ig) + &
                 fac111 * absa(ind1+10,ig)) + &
                 colh2o(iplon,lay)  * &
                 (selffac(iplon,lay)  * (selfref(inds,ig) + &
                 selffrac(iplon,lay)  * &
                 (selfref(inds+1,ig) - selfref(inds,ig))) + &
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig)))) 
            taur(iplon,lay,ngs16+ig)  = tauray
      enddo
         else


         speccomb = colh2o(iplon,lay)  + strrat*colco2(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4.*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(17) + js
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(17) + js
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng17
            taug(iplon,lay,ngs16+ig)  = speccomb * &
                (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig) + &
                 fac001 * absb(ind1,ig) + &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig)) + &
                 colh2o(iplon,lay)  * &
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
!            ssa(lay,ngs16+ig) = tauray/taug(lay,ngs16+ig)
           
            taur(iplon,lay,ngs16+ig)  = tauray
         enddo
         endif
         enddo
      enddo
!$acc end kernels
        

!$acc kernels
      do iplon = 1, ncol      
! Upper atmosphere loop
        laysolfr = nlayers 
      do lay = 2, nlayers
          if (lay > laytrop(iplon)) then 
          
        if ((jp(iplon,lay-1)  .lt. layreffr) .and. (jp(iplon,lay)  .ge. layreffr)) then
            laysolfr = lay
          end if
          
          if (lay == laysolfr) then
              
              
            speccomb = colh2o(iplon,lay)  + strrat*colco2(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1.  )
           do ig = 1, ng17 
             
                 
            if (isolvar .lt. 0) &
               sfluxzen(iplon,ngs16+ig) = sfluxref(ig,js) + &
                         fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
            if (isolvar .ge. 0 .and. isolvar .le. 2) &
               ssi(iplon,ngs16+ig) = svar_f * (facbrght(ig,js) + &
                         fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                         svar_s * (snsptdrk(ig,js) + &
                         fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                         svar_i * (irradnce(ig,js) + &
                         fs * (irradnce(ig,js+1) - irradnce(ig,js)))
            if (isolvar .eq. 3) &
               ssi(iplon,ngs16+ig) = svar_f_bnd(ngb(ngs16+ig)) * (facbrght(ig,js) + &
                         fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                         svar_s_bnd(ngb(ngs16+ig)) * (snsptdrk(ig,js) + &
                         fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                         svar_i_bnd(ngb(ngs16+ig)) * (irradnce(ig,js) + &
                         fs * (irradnce(ig,js+1) - irradnce(ig,js)))
              end do
          end if
          end if
      enddo
      enddo   
!$acc end kernels      
      end subroutine taumol17

!----------------------------------------------------------------------------
      subroutine taumol18(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng18, ngs17
      use rrsw_kg18, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray, strrat
      integer :: iplon


      strrat = 38.9589
      layreffr = 6
!$acc kernels      
      do iplon = 1, ncol
          laysolfr = laytrop(iplon)
          do lay = 1, laytrop(iplon)
              speccomb = colh2o(iplon,lay)  + strrat*colch4(iplon,lay) 
                 specparm = colh2o(iplon,lay) /speccomb 
                 if (specparm .ge. oneminus) specparm = oneminus
                specmult = 8. *(specparm)
                js = 1 + int(specmult)
                 fs = mod(specmult, 1.  )
                  if (jp(iplon,lay)  .lt. layreffr .and. jp(iplon,lay+1)  .ge. layreffr) &
            laysolfr = min(lay+1,laytrop(iplon) )
              do ig = 1, ng18
      
                   
            if (lay .eq. laysolfr .and. isolvar .lt. 0) &
               sfluxzen(iplon,ngs17+ig) = sfluxref(ig,js) + &
                         fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
            if (lay .eq. laysolfr .and. isolvar .ge. 0 .and. isolvar .le. 2) &
               ssi(iplon,ngs17+ig) = svar_f * (facbrght(ig,js) + &
                         fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                         svar_s * (snsptdrk(ig,js) + &
                         fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                         svar_i * (irradnce(ig,js) + &
                         fs * (irradnce(ig,js+1) - irradnce(ig,js)))
            if (lay .eq. laysolfr .and. isolvar .eq. 3) &
               ssi(iplon,ngs17+ig) = svar_f_bnd(ngb(ngs17+ig)) * (facbrght(ig,js) + &
                         fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                         svar_s_bnd(ngb(ngs17+ig)) * (snsptdrk(ig,js) + &
                         fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                         svar_i_bnd(ngb(ngs17+ig)) * (irradnce(ig,js) + &
                         fs * (irradnce(ig,js+1) - irradnce(ig,js)))
             end do
          end do
      end do
!$acc end kernels
      
!$acc kernels 
 do iplon = 1, ncol

      do lay = 1, nlayers
      if (lay <= laytrop(iplon)) then
          !do lay = 1, laytrop(iplon) 
       
         speccomb = colh2o(iplon,lay)  + strrat*colch4(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8.*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(18) + js
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(18) + js
         inds = indself(iplon,lay) 
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng18
            taug(iplon,lay,ngs17+ig)  = speccomb * &
                (fac000 * absa(ind0,ig) + &
                 fac100 * absa(ind0+1,ig) + &
                 fac010 * absa(ind0+9,ig) + &
                 fac110 * absa(ind0+10,ig) + &
                 fac001 * absa(ind1,ig) + &
                 fac101 * absa(ind1+1,ig) + &
                 fac011 * absa(ind1+9,ig) + &
                 fac111 * absa(ind1+10,ig)) + &
                 colh2o(iplon,lay)  * &
                 (selffac(iplon,lay)  * (selfref(inds,ig) + &
                 selffrac(iplon,lay)  * &
                 (selfref(inds+1,ig) - selfref(inds,ig))) + &
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig)))) 
!            ssa(lay,ngs17+ig) = tauray/taug(lay,ngs17+ig)
        
            taur(iplon,lay,ngs17+ig)  = tauray
      enddo

       else

! Upper atmosphere loop
              
!do lay = laytrop(iplon) +1, nlayers
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(18) + 1
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(18) + 1
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng18
            taug(iplon,lay,ngs17+ig)  = colch4(iplon,lay)  * &
                (fac00(iplon,lay)  * absb(ind0,ig) + &
                 fac10(iplon,lay)  * absb(ind0+1,ig) + &
                 fac01(iplon,lay)  * absb(ind1,ig) + &	  
                 fac11(iplon,lay)  * absb(ind1+1,ig)) 
!           ssa(lay,ngs17+ig) = tauray/taug(lay,ngs17+ig)
           taur(iplon,lay,ngs17+ig)  = tauray
         enddo
       end if
         enddo
       enddo

 
       
!$acc end kernels
       end subroutine taumol18

!----------------------------------------------------------------------------
      subroutine taumol19(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng19, ngs18
      use rrsw_kg19, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray, strrat
      integer :: iplon


	        
      strrat = 5.49281 
      layreffr = 3      
      
!$acc kernels 
 do iplon=1,ncol

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  
      laysolfr = laytrop(iplon) 
  
! Lower atmosphere loop      
      do lay = 1, laytrop(iplon) 
            
        if (jp(iplon,lay)  .lt. layreffr .and. jp(iplon,lay+1)  .ge. layreffr) &
            laysolfr = min(lay+1,laytrop(iplon) )
     
         if (lay .eq. laysolfr) then 
                 speccomb = colh2o(iplon,lay)  + strrat*colco2(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1.  )
        
         do ig = 1 , ng19
            if (isolvar .lt. 0) &
               sfluxzen(iplon,ngs18+ig) = sfluxref(ig,js) + &
                         fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
            if (isolvar .ge. 0 .and. isolvar .le. 2) &
               ssi(iplon,ngs18+ig) = svar_f * (facbrght(ig,js) + &
                         fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                         svar_s * (snsptdrk(ig,js) + &
                         fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                         svar_i * (irradnce(ig,js) + &
                         fs * (irradnce(ig,js+1) - irradnce(ig,js)))
            if (isolvar .eq. 3) &
               ssi(iplon,ngs18+ig) = svar_f_bnd(ngb(ngs18+ig)) * (facbrght(ig,js) + &
                         fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                         svar_s_bnd(ngb(ngs18+ig)) * (snsptdrk(ig,js) + &
                         fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                         svar_i_bnd(ngb(ngs18+ig)) * (irradnce(ig,js) + &
                         fs * (irradnce(ig,js+1) - irradnce(ig,js)))
         end do
         end if
      end do
 end do
!$acc end kernels
      
      

      
!$acc kernels 
 do iplon=1,ncol

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  



! Lower atmosphere loop      
      do lay = 1, nlayers
          if (lay <= laytrop(iplon)) then
       
         speccomb = colh2o(iplon,lay)  + strrat*colco2(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8.*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(19) + js
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(19) + js
         inds = indself(iplon,lay) 
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1 , ng19
            taug(iplon,lay,ngs18+ig)  = speccomb * &
                (fac000 * absa(ind0,ig) + &
                 fac100 * absa(ind0+1,ig) + &
                 fac010 * absa(ind0+9,ig) + &
                 fac110 * absa(ind0+10,ig) + &
                 fac001 * absa(ind1,ig) + &
                 fac101 * absa(ind1+1,ig) + &
                 fac011 * absa(ind1+9,ig) + &
                 fac111 * absa(ind1+10,ig)) + &
                 colh2o(iplon,lay)  * &
                 (selffac(iplon,lay)  * (selfref(inds,ig) + &
                 selffrac(iplon,lay)  * &
                 (selfref(inds+1,ig) - selfref(inds,ig))) + & 
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig)))) 
!            ssa(lay,ngs18+ig) = tauray/taug(lay,ngs18+ig)
            taur(iplon,lay,ngs18+ig)  = tauray   
      enddo
     else

! Upper atmosphere loop
  
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(19) + 1
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(19) + 1
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1 , ng19
            taug(iplon,lay,ngs18+ig)  = colco2(iplon,lay)  * &
                (fac00(iplon,lay)  * absb(ind0,ig) + &
                 fac10(iplon,lay)  * absb(ind0+1,ig) + &
                 fac01(iplon,lay)  * absb(ind1,ig) + &
                 fac11(iplon,lay)  * absb(ind1+1,ig)) 
!            ssa(lay,ngs18+ig) = tauray/taug(lay,ngs18+ig) 
            taur(iplon,lay,ngs18+ig)  = tauray   
         enddo
         end if
         enddo
      enddo

      
!$acc end kernels
      end subroutine taumol19

!----------------------------------------------------------------------------
      subroutine taumol20(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng20, ngs19
      use rrsw_kg20, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, absch4, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb


      implicit none

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray
        integer :: iplon

      layreffr = 3
!$acc kernels loop independent private(laysolfr)
do iplon = 1, ncol
    laysolfr = laytrop(iplon)
    do lay = 1, laytrop(iplon)
        if (jp(iplon,lay)  .lt. layreffr .and. jp(iplon,lay+1)  .ge. layreffr) &
            laysolfr = min(lay+1,laytrop(iplon) )
         if (lay .eq. laysolfr) then 
             do ig = 1, ng20 
                 if (isolvar .lt. 0) &
                     sfluxzen(iplon,ngs19+ig) = sfluxref(ig) 
                 if (isolvar .ge. 0 .and. isolvar .le. 2) &
                     ssi(iplon,ngs19+ig) = svar_f * facbrght(ig) + &
                              svar_s * snsptdrk(ig) + &
                              svar_i * irradnce(ig)
                 if (isolvar .eq. 3) &
                     ssi(iplon,ngs19+ig) = svar_f_bnd(ngb(ngs19+ig)) * facbrght(ig) + &
                              svar_s_bnd(ngb(ngs19+ig)) * snsptdrk(ig) + &
                              svar_i_bnd(ngb(ngs19+ig)) * irradnce(ig)
             end do
         end if
    end do
end do
!$acc end kernels

!$acc kernels 
 do iplon=1,ncol

      do lay = 1, nlayers 
          if (lay <= laytrop(iplon)) then
         
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(20) + 1
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(20) + 1
         inds = indself(iplon,lay) 
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng20
            taug(iplon,lay,ngs19+ig)  = colh2o(iplon,lay)  * &
               ((fac00(iplon,lay)  * absa(ind0,ig) + &
                 fac10(iplon,lay)  * absa(ind0+1,ig) + &
                 fac01(iplon,lay)  * absa(ind1,ig) + &
                 fac11(iplon,lay)  * absa(ind1+1,ig)) + &
                 selffac(iplon,lay)  * (selfref(inds,ig) + & 
                 selffrac(iplon,lay)  * &
                 (selfref(inds+1,ig) - selfref(inds,ig))) + &
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig)))) &
                 + colch4(iplon,lay)  * absch4(ig)
!            ssa(lay,ngs19+ig) = tauray/taug(lay,ngs19+ig)
            taur(iplon,lay,ngs19+ig)  = tauray 
           
      enddo
      else

! Upper atmosphere loop
      
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(20) + 1
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(20) + 1
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng20
            taug(iplon,lay,ngs19+ig)  = colh2o(iplon,lay)  * &
                (fac00(iplon,lay)  * absb(ind0,ig) + &
                 fac10(iplon,lay)  * absb(ind0+1,ig) + &
                 fac01(iplon,lay)  * absb(ind1,ig) + &
                 fac11(iplon,lay)  * absb(ind1+1,ig) + &
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig)))) + &
                 colch4(iplon,lay)  * absch4(ig)
!            ssa(lay,ngs19+ig) = tauray/taug(lay,ngs19+ig)
            taur(iplon,lay,ngs19+ig)  = tauray 
         enddo
         end if
         enddo
      enddo

      
!$acc end kernels
      end subroutine taumol20

!----------------------------------------------------------------------------
      subroutine taumol21(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng21, ngs20
      use rrsw_kg21, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray, strrat
        integer :: iplon
        

      strrat = 0.0045321 
      layreffr = 8
     
        
!$acc kernels loop independent private(laysolfr)
do iplon=1,ncol
     laysolfr = laytrop(iplon)        
!$acc loop seq
    do lay=1,laytrop(iplon)
         if (jp(iplon,lay)  .lt. layreffr .and. jp(iplon,lay+1)  .ge. layreffr) &
            laysolfr = min(lay+1,laytrop(iplon) )
          if (lay .eq. laysolfr) then 
                speccomb = colh2o(iplon,lay)  + strrat*colco2(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1.  )
               do ig = 1, ng21
               if (isolvar .lt. 0) &
                  sfluxzen(iplon,ngs20+ig) = sfluxref(ig,js) + &
                           fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
               if (isolvar .ge. 0 .and. isolvar .le. 2) &
                  ssi(iplon,ngs20+ig) = svar_f * (facbrght(ig,js) + &
                           fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                           svar_s * (snsptdrk(ig,js) + &
                           fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                           svar_i * (irradnce(ig,js) + &
                           fs * (irradnce(ig,js+1) - irradnce(ig,js)))
               if (isolvar .eq. 3) &
                  ssi(iplon,ngs20+ig) = svar_f_bnd(ngb(ngs20+ig)) * (facbrght(ig,js) + &
                           fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                           svar_s_bnd(ngb(ngs20+ig)) * (snsptdrk(ig,js) + &
                           fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                           svar_i_bnd(ngb(ngs20+ig)) * (irradnce(ig,js) + &
                           fs * (irradnce(ig,js+1) - irradnce(ig,js)))
               end do
          end if
    end do
end do        
!$acc end kernels

       
!$acc kernels 
 do iplon=1,ncol  
! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  

      
! Lower atmosphere loop
      do lay = 1, nlayers 
      if (lay <= laytrop(iplon)) then 
        
         speccomb = colh2o(iplon,lay)  + strrat*colco2(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8.*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(21) + js
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(21) + js
         inds = indself(iplon,lay) 
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng21
            taug(iplon,lay,ngs20+ig)  = speccomb * &
                (fac000 * absa(ind0,ig) + &
                 fac100 * absa(ind0+1,ig) + &
                 fac010 * absa(ind0+9,ig) + &
                 fac110 * absa(ind0+10,ig) + &
                 fac001 * absa(ind1,ig) + &
                 fac101 * absa(ind1+1,ig) + &
                 fac011 * absa(ind1+9,ig) + &
                 fac111 * absa(ind1+10,ig)) + &
                 colh2o(iplon,lay)  * &
                 (selffac(iplon,lay)  * (selfref(inds,ig) + &
                 selffrac(iplon,lay)  * &
                 (selfref(inds+1,ig) - selfref(inds,ig))) + &
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig))))
!            ssa(lay,ngs20+ig) = tauray/taug(lay,ngs20+ig)
          
            taur(iplon,lay,ngs20+ig)  = tauray
      enddo
      else

! Upper atmosphere loop

         speccomb = colh2o(iplon,lay)  + strrat*colco2(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4.*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(21) + js
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(21) + js
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng21
            taug(iplon,lay,ngs20+ig)  = speccomb * &
                (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig) + &
                 fac001 * absb(ind1,ig) + &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig)) + &
                 colh2o(iplon,lay)  * &
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig)))
!            ssa(lay,ngs20+ig) = tauray/taug(lay,ngs20+ig)
            taur(iplon,lay,ngs20+ig)  = tauray
         enddo
         end if
         enddo
      enddo

      
!$acc end kernels
      end subroutine taumol21

!----------------------------------------------------------------------------
      subroutine taumol22(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng22, ngs21
      use rrsw_kg22, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray, o2adj, o2cont, strrat
        integer :: iplon

! The following factor is the ratio of total O2 band intensity (lines 
! and Mate continuum) to O2 band intensity (line only).  It is needed
! to adjust the optical depths since the k's include only lines.
      o2adj = 1.6
      
! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  

      strrat = 0.022708
      layreffr = 2
      
!$acc kernels loop independent private(laysolfr)
 do iplon=1,ncol

      laysolfr = laytrop(iplon) 

! Lower atmosphere loop
!$acc loop seq
      do lay = 1, laytrop(iplon) 
                 if (jp(iplon,lay)  .lt. layreffr .and. jp(iplon,lay+1)  .ge. layreffr) &
            laysolfr = min(lay+1,laytrop(iplon) )
                 
                if (lay .eq. laysolfr) then 
            speccomb = colh2o(iplon,lay)  + o2adj*strrat*colo2(iplon,lay) 
            specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8.*(specparm)
!         odadj = specparm + o2adj * (1. - specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
                do ig = 1, ng22                                 
                                 
                if (isolvar .lt. 0) &
                     sfluxzen(iplon,ngs21+ig) = sfluxref(ig,js) + &
                              fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
                if (isolvar .ge. 0 .and. isolvar .le. 2) &
                     ssi(iplon,ngs21+ig) = svar_f * (facbrght(ig,js) + &
                              fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                              svar_s * (snsptdrk(ig,js) + &
                              fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                              svar_i * (irradnce(ig,js) + &
                              fs * (irradnce(ig,js+1) - irradnce(ig,js)))
                if (isolvar .eq. 3) &
                     ssi(iplon,ngs21+ig) = svar_f_bnd(ngb(ngs21+ig)) * (facbrght(ig,js) + &
                              fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                              svar_s_bnd(ngb(ngs21+ig)) * (snsptdrk(ig,js) + &
                              fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                              svar_i_bnd(ngb(ngs21+ig)) * (irradnce(ig,js) + &
                              fs * (irradnce(ig,js+1) - irradnce(ig,js)))
                end do
                end if
      end do
 end do
 !$acc end kernels
 
!$acc kernels 
 do iplon=1,ncol

      laysolfr = laytrop(iplon) 

! Lower atmosphere loop
      do lay = 1, nlayers 
          if (lay<=laytrop(iplon)) then
  
         o2cont = 4.35e-4 *colo2(iplon,lay) /(350.0 *2.0 )
         speccomb = colh2o(iplon,lay)  + o2adj*strrat*colo2(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
!         odadj = specparm + o2adj * (1.  - specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1.  )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(22) + js
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(22) + js
         inds = indself(iplon,lay) 
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng22
            taug(iplon,lay,ngs21+ig)  = speccomb * &
                (fac000 * absa(ind0,ig) + &
                 fac100 * absa(ind0+1,ig) + &
                 fac010 * absa(ind0+9,ig) + &
                 fac110 * absa(ind0+10,ig) + &
                 fac001 * absa(ind1,ig) + &
                 fac101 * absa(ind1+1,ig) + &
                 fac011 * absa(ind1+9,ig) + &
                 fac111 * absa(ind1+10,ig)) + &
                 colh2o(iplon,lay)  * &
                 (selffac(iplon,lay)  * (selfref(inds,ig) + &
                 selffrac(iplon,lay)  * &
                  (selfref(inds+1,ig) - selfref(inds,ig))) + &
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig)))) &
                 + o2cont
!            ssa(lay,ngs21+ig) = tauray/taug(lay,ngs21+ig)

            taur(iplon,lay,ngs21+ig)  = tauray
      enddo
      else

! Upper atmosphere loop
      
         o2cont = 4.35e-4 *colo2(iplon,lay) /(350.0 *2.0 )
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(22) + 1
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(22) + 1
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng22
            taug(iplon,lay,ngs21+ig)  = colo2(iplon,lay)  * o2adj * &
                (fac00(iplon,lay)  * absb(ind0,ig) + &
                 fac10(iplon,lay)  * absb(ind0+1,ig) + &
                 fac01(iplon,lay)  * absb(ind1,ig) + &
                 fac11(iplon,lay)  * absb(ind1+1,ig)) + &
                 o2cont
!            ssa(lay,ngs21+ig) = tauray/taug(lay,ngs21+ig)
            taur(iplon,lay,ngs21+ig)  = tauray
         enddo
         end if
         enddo
      enddo

!$acc end kernels
      end subroutine taumol22

!----------------------------------------------------------------------------
      subroutine taumol23(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng23, ngs22
      use rrsw_kg23, only : absa, ka, forref, selfref, &
                            sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray, givfac
        integer :: iplon


! Average Giver et al. correction factor for this band.
      givfac = 1.029

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  

      layreffr = 6
      
!$acc kernels loop independent private(laysolfr)
 do iplon=1,ncol

      laysolfr = laytrop(iplon) 

! Lower atmosphere loop
!$acc loop seq
      do lay = 1, laytrop(iplon) 
         if (jp(iplon,lay)  .lt. layreffr .and. jp(iplon,lay+1)  .ge. layreffr) &
            laysolfr = min(lay+1,laytrop(iplon) )

          if (lay .eq. laysolfr) then 
         do ig = 1, ng23
            if (isolvar .lt. 0) &
               sfluxzen(iplon,ngs22+ig) = sfluxref(ig) 
            if (isolvar .ge. 0 .and. isolvar .le. 2) &
               ssi(iplon,ngs22+ig) = svar_f * facbrght(ig) + &
                         svar_s * snsptdrk(ig) + &
                         svar_i * irradnce(ig)
            if (isolvar .eq. 3) &
               ssi(iplon,ngs22+ig) = svar_f_bnd(ngb(ngs22+ig)) * facbrght(ig) + &
                         svar_s_bnd(ngb(ngs22+ig)) * snsptdrk(ig) + &
                         svar_i_bnd(ngb(ngs22+ig)) * irradnce(ig)
            end do
          end if
      end do
end do      
!$acc end kernels   
      
      
!$acc kernels 
 do iplon=1,ncol

     

! Lower atmosphere loop
      do lay = 1, nlayers 
         if (lay <= laytrop(iplon)) then
         if (jp(iplon,lay)  .lt. layreffr .and. jp(iplon,lay+1)  .ge. layreffr) &
            laysolfr = min(lay+1,laytrop(iplon) )
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(23) + 1
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(23) + 1
         inds = indself(iplon,lay) 
         indf = indfor(iplon,lay) 

         do ig = 1, ng23
            tauray = colmol(iplon,lay)  * rayl(ig)
            taug(iplon,lay,ngs22+ig)  = colh2o(iplon,lay)  * &
                (givfac * (fac00(iplon,lay)  * absa(ind0,ig) + &
                 fac10(iplon,lay)  * absa(ind0+1,ig) + &
                 fac01(iplon,lay)  * absa(ind1,ig) + &
                 fac11(iplon,lay)  * absa(ind1+1,ig)) + &
                 selffac(iplon,lay)  * (selfref(inds,ig) + &
                 selffrac(iplon,lay)  * &
                 (selfref(inds+1,ig) - selfref(inds,ig))) + &
                 forfac(iplon,lay)  * (forref(indf,ig) + &
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig)))) 
!            ssa(lay,ngs22+ig) = tauray/taug(lay,ngs22+ig)
           
            taur(iplon,lay,ngs22+ig)  = tauray
      enddo
      else

! Upper atmosphere loop
      
         do ig = 1, ng23
!            taug(lay,ngs22+ig) = colmol(lay) * rayl(ig)
!            ssa(lay,ngs22+ig) = 1.0
            taug(iplon,lay,ngs22+ig)  = 0. 
            taur(iplon,lay,ngs22+ig)  = colmol(iplon,lay)  * rayl(ig) 
         enddo
      end if
      enddo
      enddo


!$acc end kernels
      end subroutine taumol23

!----------------------------------------------------------------------------
      subroutine taumol24(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng24, ngs23
      use rrsw_kg24, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, abso3a, abso3b, rayla, raylb, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray, strrat
        integer :: iplon

        
      strrat = 0.124692 
      layreffr = 1   
        
!$acc kernels loop independent private(laysolfr)
 do iplon=1,ncol
! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  


      laysolfr = laytrop(iplon) 

! Lower atmosphere loop
!$acc loop independent
      do lay = 1, laytrop(iplon) 
                   if (jp(iplon,lay)  .lt. layreffr .and. jp(iplon,lay+1)  .ge. layreffr) &
            laysolfr = min(lay+1,laytrop(iplon) )
          if (lay .eq. laysolfr) then
                 speccomb = colh2o(iplon,lay)  + strrat*colo2(iplon,lay) 
            specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8.*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
          do ig = 1, ng24
            if (isolvar .lt. 0) &
               sfluxzen(iplon,ngs23+ig) = sfluxref(ig,js) + &
                         fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
            if (isolvar .ge. 0 .and. isolvar .le. 2) &
               ssi(iplon,ngs23+ig) = svar_f * (facbrght(ig,js) + &
                         fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                         svar_s * (snsptdrk(ig,js) + &
                         fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                         svar_i * (irradnce(ig,js) + &
                         fs * (irradnce(ig,js+1) - irradnce(ig,js)))
            if (isolvar .eq. 3) &
               ssi(iplon,ngs23+ig) = svar_f_bnd(ngb(ngs23+ig)) * (facbrght(ig,js) + &
                         fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                         svar_s_bnd(ngb(ngs23+ig)) * (snsptdrk(ig,js) + &
                         fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                         svar_i_bnd(ngb(ngs23+ig)) * (irradnce(ig,js) + &
                         fs * (irradnce(ig,js+1) - irradnce(ig,js)))
          end do
          end if
      end do
 end do
 !$acc end kernels
     
        
!$acc kernels 
 do iplon=1,ncol
! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  




! Lower atmosphere loop
      do lay = 1, nlayers 
          if (lay <= laytrop(iplon)) then

         speccomb = colh2o(iplon,lay)  + strrat*colo2(iplon,lay) 
         specparm = colh2o(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1.  )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(24) + js
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(24) + js
         inds = indself(iplon,lay) 
         indf = indfor(iplon,lay) 

         do ig = 1, ng24
            tauray = colmol(iplon,lay)  * (rayla(ig,js) + &
               fs * (rayla(ig,js+1) - rayla(ig,js)))
            taug(iplon,lay,ngs23+ig)  = speccomb * &
                (fac000 * absa(ind0,ig) + &
                 fac100 * absa(ind0+1,ig) + &
                 fac010 * absa(ind0+9,ig) + &
                 fac110 * absa(ind0+10,ig) + &
                 fac001 * absa(ind1,ig) + &
                 fac101 * absa(ind1+1,ig) + &
                 fac011 * absa(ind1+9,ig) + &
                 fac111 * absa(ind1+10,ig)) + &
                 colo3(iplon,lay)  * abso3a(ig) + &
                 colh2o(iplon,lay)  * & 
                 (selffac(iplon,lay)  * (selfref(inds,ig) + &
                 selffrac(iplon,lay)  * &
                 (selfref(inds+1,ig) - selfref(inds,ig))) + &
                 forfac(iplon,lay)  * (forref(indf,ig) + & 
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig))))
!            ssa(lay,ngs23+ig) = tauray/taug(lay,ngs23+ig)
           
            taur(iplon,lay,ngs23+ig)  = tauray
      enddo
     else

! Upper atmosphere loop
      
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(24) + 1
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(24) + 1

         do ig = 1, ng24
            tauray = colmol(iplon,lay)  * raylb(ig)
            taug(iplon,lay,ngs23+ig)  = colo2(iplon,lay)  * &
                (fac00(iplon,lay)  * absb(ind0,ig) + &
                 fac10(iplon,lay)  * absb(ind0+1,ig) + &
                 fac01(iplon,lay)  * absb(ind1,ig) + &
                 fac11(iplon,lay)  * absb(ind1+1,ig)) + &
                 colo3(iplon,lay)  * abso3b(ig)
!            ssa(lay,ngs23+ig) = tauray/taug(lay,ngs23+ig)
            taur(iplon,lay,ngs23+ig)  = tauray
         enddo
         endif
         enddo
      enddo

!$acc end kernels
      end subroutine taumol24

!----------------------------------------------------------------------------
      subroutine taumol25(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng25, ngs24
      use rrsw_kg25, only : absa, ka, &
                            sfluxref, abso3a, abso3b, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray
    integer :: iplon

       
!$acc kernels 
 do iplon=1,ncol
! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  

      layreffr = 2
      laysolfr = laytrop(iplon) 

! Lower atmosphere loop
      do lay = 1, laytrop(iplon) 
         if (jp(iplon,lay)  .lt. layreffr .and. jp(iplon,lay+1)  .ge. layreffr) &
            laysolfr = min(lay+1,laytrop(iplon) )
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(25) + 1
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(25) + 1

         do ig = 1, ng25
            tauray = colmol(iplon,lay)  * rayl(ig)
            taug(iplon,lay,ngs24+ig)  = colh2o(iplon,lay)  * &
                (fac00(iplon,lay)  * absa(ind0,ig) + &
                 fac10(iplon,lay)  * absa(ind0+1,ig) + &
                 fac01(iplon,lay)  * absa(ind1,ig) + &
                 fac11(iplon,lay)  * absa(ind1+1,ig)) + &
                 colo3(iplon,lay)  * abso3a(ig) 
!            ssa(lay,ngs24+ig) = tauray/taug(lay,ngs24+ig)
            if (lay .eq. laysolfr .and. isolvar .lt. 0) &
               sfluxzen(iplon,ngs24+ig) = sfluxref(ig) 
            if (lay .eq. laysolfr .and. isolvar .ge. 0 .and. isolvar .le. 2) &
               ssi(iplon,ngs24+ig) = svar_f * facbrght(ig) + &
                         svar_s * snsptdrk(ig) + &
                         svar_i * irradnce(ig)
            if (lay .eq. laysolfr .and. isolvar .eq. 3) &
               ssi(iplon,ngs24+ig) = svar_f_bnd(ngb(ngs24+ig)) * facbrght(ig) + &
                         svar_s_bnd(ngb(ngs24+ig)) * snsptdrk(ig) + &
                         svar_i_bnd(ngb(ngs24+ig)) * irradnce(ig)
            taur(iplon,lay,ngs24+ig)  = tauray
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop(iplon) +1, nlayers
         do ig = 1, ng25
            tauray = colmol(iplon,lay)  * rayl(ig)
            taug(iplon,lay,ngs24+ig)  = colo3(iplon,lay)  * abso3b(ig) 
!            ssa(lay,ngs24+ig) = tauray/taug(lay,ngs24+ig)
            taur(iplon,lay,ngs24+ig)  = tauray
         enddo
         enddo
      enddo

!$acc end kernels
      end subroutine taumol25

!----------------------------------------------------------------------------
      subroutine taumol26(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng26, ngs25
      use rrsw_kg26, only : sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray
    integer :: iplon

       
!$acc kernels 
 do iplon=1,ncol
! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  

      laysolfr = laytrop(iplon) 

! Lower atmosphere loop
      do lay = 1, laytrop(iplon) 
         do ig = 1, ng26 
!            taug(lay,ngs25+ig) = colmol(lay) * rayl(ig)
!            ssa(lay,ngs25+ig) = 1.0
            if (lay .eq. laysolfr .and. isolvar .lt. 0) &
               sfluxzen(iplon,ngs25+ig) = sfluxref(ig) 
            if (lay .eq. laysolfr .and. isolvar .ge. 0 .and. isolvar .le. 2) &
               ssi(iplon,ngs25+ig) = svar_f * facbrght(ig) + &
                         svar_s * snsptdrk(ig) + &
                         svar_i * irradnce(ig)
            if (lay .eq. laysolfr .and. isolvar .eq. 3) &
               ssi(iplon,ngs25+ig) = svar_f_bnd(ngb(ngs25+ig)) * facbrght(ig) + &
                         svar_s_bnd(ngb(ngs25+ig)) * snsptdrk(ig) + &
                         svar_i_bnd(ngb(ngs25+ig)) * irradnce(ig)
            taug(iplon,lay,ngs25+ig)  = 0. 
            taur(iplon,lay,ngs25+ig)  = colmol(iplon,lay)  * rayl(ig) 
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop(iplon) +1, nlayers
         do ig = 1, ng26
!            taug(lay,ngs25+ig) = colmol(lay) * rayl(ig)
!            ssa(lay,ngs25+ig) = 1.0
            taug(iplon,lay,ngs25+ig)  = 0. 
            taur(iplon,lay,ngs25+ig)  = colmol(iplon,lay)  * rayl(ig) 
         enddo
         enddo
      enddo

!$acc end kernels
      end subroutine taumol26

!----------------------------------------------------------------------------
      subroutine taumol27(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng27, ngs26
      use rrsw_kg27, only : absa, ka, absb, kb, sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray, scalekur
    integer :: iplon


!$acc kernels 
 do iplon=1,ncol
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
      do lay = 1, laytrop(iplon) 
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(27) + 1
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(27) + 1

         do ig = 1, ng27
            tauray = colmol(iplon,lay)  * rayl(ig)
            taug(iplon,lay,ngs26+ig)  = colo3(iplon,lay)  * &
                (fac00(iplon,lay)  * absa(ind0,ig) + &
                 fac10(iplon,lay)  * absa(ind0+1,ig) + &
                 fac01(iplon,lay)  * absa(ind1,ig) + &
                 fac11(iplon,lay)  * absa(ind1+1,ig))
!            ssa(lay,ngs26+ig) = tauray/taug(lay,ngs26+ig)
            taur(iplon,lay,ngs26+ig)  = tauray
         enddo
      enddo

      laysolfr = nlayers

! Upper atmosphere loop
      do lay = laytrop(iplon) +1, nlayers
         if (jp(iplon,lay-1)  .lt. layreffr .and. jp(iplon,lay)  .ge. layreffr) &
            laysolfr = lay
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(27) + 1
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(27) + 1

         do ig = 1, ng27
            tauray = colmol(iplon,lay)  * rayl(ig)
            taug(iplon,lay,ngs26+ig)  = colo3(iplon,lay)  * &
                (fac00(iplon,lay)  * absb(ind0,ig) + &
                 fac10(iplon,lay)  * absb(ind0+1,ig) + &
                 fac01(iplon,lay)  * absb(ind1,ig) + & 
                 fac11(iplon,lay)  * absb(ind1+1,ig))
!            ssa(lay,ngs26+ig) = tauray/taug(lay,ngs26+ig)
            if (lay .eq. laysolfr .and. isolvar .lt. 0) &
               sfluxzen(iplon,ngs26+ig) = scalekur * sfluxref(ig) 
            if (lay .eq. laysolfr .and. isolvar .ge. 0 .and. isolvar .le. 2) &
               ssi(iplon,ngs26+ig) = svar_f * facbrght(ig) + &
                         svar_s * snsptdrk(ig) + &
                         svar_i * irradnce(ig)
            if (lay .eq. laysolfr .and. isolvar .eq. 3) &
               ssi(iplon,ngs26+ig) = svar_f_bnd(ngb(ngs26+ig)) * facbrght(ig) + &
                         svar_s_bnd(ngb(ngs26+ig)) * snsptdrk(ig) + &
                         svar_i_bnd(ngb(ngs26+ig)) * irradnce(ig)
            taur(iplon,lay,ngs26+ig)  = tauray
         enddo
         enddo
      enddo

!$acc end kernels
      end subroutine taumol27

!----------------------------------------------------------------------------
      subroutine taumol28(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng28, ngs27
      use rrsw_kg28, only : absa, ka, absb, kb, sfluxref, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray, strrat
        integer :: iplon

       
!$acc kernels 
 do iplon=1,ncol
! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  

      strrat = 6.67029e-07
      layreffr = 58

! Lower atmosphere loop
      do lay = 1, laytrop(iplon) 
         speccomb = colo3(iplon,lay)  + strrat*colo2(iplon,lay) 
         specparm = colo3(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8.*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(28) + js
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(28) + js
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng28
            taug(iplon,lay,ngs27+ig)  = speccomb * &
                (fac000 * absa(ind0,ig) + &
                 fac100 * absa(ind0+1,ig) + &
                 fac010 * absa(ind0+9,ig) + &
                 fac110 * absa(ind0+10,ig) + &
                 fac001 * absa(ind1,ig) + &
                 fac101 * absa(ind1+1,ig) + &
                 fac011 * absa(ind1+9,ig) + &
                 fac111 * absa(ind1+10,ig)) 
!            ssa(lay,ngs27+ig) = tauray/taug(lay,ngs27+ig)
            taur(iplon,lay,ngs27+ig)  = tauray
         enddo
      enddo

      laysolfr = nlayers

! Upper atmosphere loop
      do lay = laytrop(iplon) +1, nlayers
         if (jp(iplon,lay-1)  .lt. layreffr .and. jp(iplon,lay)  .ge. layreffr) &
            laysolfr = lay
         speccomb = colo3(iplon,lay)  + strrat*colo2(iplon,lay) 
         specparm = colo3(iplon,lay) /speccomb 
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4.*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult, 1. )
         fac000 = (1.  - fs) * fac00(iplon,lay) 
         fac010 = (1.  - fs) * fac10(iplon,lay) 
         fac100 = fs * fac00(iplon,lay) 
         fac110 = fs * fac10(iplon,lay) 
         fac001 = (1.  - fs) * fac01(iplon,lay) 
         fac011 = (1.  - fs) * fac11(iplon,lay) 
         fac101 = fs * fac01(iplon,lay) 
         fac111 = fs * fac11(iplon,lay) 
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(28) + js
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(28) + js
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng28
            taug(iplon,lay,ngs27+ig)  = speccomb * &
                (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig) + &
                 fac001 * absb(ind1,ig) + &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig)) 
!            ssa(lay,ngs27+ig) = tauray/taug(lay,ngs27+ig)
            if (lay .eq. laysolfr .and. isolvar .lt. 0) &
               sfluxzen(iplon,ngs27+ig) = sfluxref(ig,js) + &
                         fs * (sfluxref(ig,js+1) - sfluxref(ig,js))
            if (lay .eq. laysolfr .and. isolvar .ge. 0 .and. isolvar .le. 2) &
               ssi(iplon,ngs27+ig) = svar_f * (facbrght(ig,js) + &
                         fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                         svar_s * (snsptdrk(ig,js) + &
                         fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                         svar_i * (irradnce(ig,js) + &
                         fs * (irradnce(ig,js+1) - irradnce(ig,js)))
            if (lay .eq. laysolfr .and. isolvar .eq. 3) &
               ssi(iplon,ngs27+ig) = svar_f_bnd(ngb(ngs27+ig)) * (facbrght(ig,js) + &
                         fs * (facbrght(ig,js+1) - facbrght(ig,js))) + &
                         svar_s_bnd(ngb(ngs27+ig)) * (snsptdrk(ig,js) + &
                         fs * (snsptdrk(ig,js+1) - snsptdrk(ig,js))) + &
                         svar_i_bnd(ngb(ngs27+ig)) * (irradnce(ig,js) + &
                         fs * (irradnce(ig,js+1) - irradnce(ig,js)))
            taur(iplon,lay,ngs27+ig)  = tauray
         enddo
         enddo
      enddo

!$acc end kernels
      end subroutine taumol28

!----------------------------------------------------------------------------
      subroutine taumol29(ncol, nlayers, &
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

! ------- Modules -------

      use parrrsw, only : ng29, ngs28
      use rrsw_kg29, only : absa, ka, absb, kb, forref, selfref, &
                            sfluxref, absh2o, absco2, rayl, &
                            irradnce, facbrght, snsptdrk
      use rrsw_wvn, only : ngb

! ------- Declarations -------
 integer , intent(in) :: ncol
      integer , intent(in) :: nlayers            ! total number of layers

      integer , intent(in) :: laytrop(:)            ! tropopause layer index
      integer , intent(in) :: jp(:,:)               ! 
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:)               !
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:)              !
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: colh2o(:,:)              ! column amount (h2o)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:)              ! column amount (co2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:)               ! column amount (o3)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:)              ! column amount (ch4)
                                                         !   Dimensions: (nlayers)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:)               ! column amount (o2)
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:)              ! 
                                                         !   Dimensions: (nlayers)

      integer , intent(in) :: indself(:,:)     
                                                         !   Dimensions: (nlayers)
      integer , intent(in) :: indfor(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfac(:,:) 
                                                         !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                         !   Dimensions: (nlayers)

      real , intent(in) :: &                     !
                       fac00(:,:) , fac01(:,:) , &             !   Dimensions: (nlayers)
                       fac10(:,:) , fac11(:,:)  

! Solar variability
      integer, intent(in) :: isolvar            ! Flag for solar variability method
      real, intent(in) :: svar_f                ! Solar variability facular multiplier
      real, intent(in) :: svar_s                ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(:)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(:)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(:)         ! Solar variability baseline irradiance multiplier (by band)

! ----- Output -----
      real, intent(out) gpu_device :: ssi(:,:)                ! spectral solar intensity with solar variability
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: sfluxzen(:,:)           ! solar source function
                                                         !   Dimensions: (ngptsw)
      real, intent(out) gpu_device :: taug(:,:,:)             ! gaseous optical depth 
                                                         !   Dimensions: (nlayers,ngptsw)
      real, intent(out) gpu_device :: taur(:,:,:)             ! Rayleigh 
                                                         !   Dimensions: (nlayers,ngptsw)
! Local

      integer :: ig, ind0, ind1, inds, indf, js, lay, laysolfr, &
                          layreffr
      real :: fac000, fac001, fac010, fac011, fac100, fac101, &
                       fac110, fac111, fs, speccomb, specmult, specparm, &
                       tauray
        integer :: iplon

            layreffr = 49  
        
!$acc kernels loop independent private (laysolfr)
 do iplon=1,ncol
        
        laysolfr = nlayers
!$acc loop seq
        do lay = laytrop(iplon) +1, nlayers
         if (jp(iplon,lay-1)  .lt. layreffr .and. jp(iplon,lay)  .ge. layreffr) &
            laysolfr = lay

            if (lay .eq. laysolfr) then 
                do ig = 1, ng29
                if (isolvar .lt. 0) &
                     sfluxzen(iplon,ngs28+ig) = sfluxref(ig) 
                if (isolvar .ge. 0 .and. isolvar .le. 2) &
                     ssi(iplon,ngs28+ig) = svar_f * facbrght(ig) + &
                              svar_s * snsptdrk(ig) + &
                              svar_i * irradnce(ig)
                if (isolvar .eq. 3) &
                     ssi(iplon,ngs28+ig) = svar_f_bnd(ngb(ngs28+ig)) * facbrght(ig) + &
                              svar_s_bnd(ngb(ngs28+ig)) * snsptdrk(ig) + &
                              svar_i_bnd(ngb(ngs28+ig)) * irradnce(ig)
                end do
            end if
        end do
 end do
 !$acc end kernels

!$acc kernels 
 do iplon=1,ncol
! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below LAYTROP, the water
! vapor self-continuum is interpolated (in temperature) separately.  


! Lower atmosphere loop
      do lay = 1, nlayers 
          if (lay <= laytrop(iplon)) then
         ind0 = ((jp(iplon,lay) -1)*5+(jt(iplon,lay) -1))*nspa(29) + 1
         ind1 = (jp(iplon,lay) *5+(jt1(iplon,lay) -1))*nspa(29) + 1
         inds = indself(iplon,lay) 
         indf = indfor(iplon,lay) 
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng29
            taug(iplon,lay,ngs28+ig)  = colh2o(iplon,lay)  * &
               ((fac00(iplon,lay)  * absa(ind0,ig) + &
                 fac10(iplon,lay)  * absa(ind0+1,ig) + &
                 fac01(iplon,lay)  * absa(ind1,ig) + &
                 fac11(iplon,lay)  * absa(ind1+1,ig)) + &
                 selffac(iplon,lay)  * (selfref(inds,ig) + &
                 selffrac(iplon,lay)  * &
                 (selfref(inds+1,ig) - selfref(inds,ig))) + &
                 forfac(iplon,lay)  * (forref(indf,ig) + & 
                 forfrac(iplon,lay)  * &
                 (forref(indf+1,ig) - forref(indf,ig)))) &
                 + colco2(iplon,lay)  * absco2(ig) 
!            ssa(lay,ngs28+ig) = tauray/taug(lay,ngs28+ig)
            taur(iplon,lay,ngs28+ig)  = tauray
      enddo
    else 



! Upper atmosphere loop
      
        
         ind0 = ((jp(iplon,lay) -13)*5+(jt(iplon,lay) -1))*nspb(29) + 1
         ind1 = ((jp(iplon,lay) -12)*5+(jt1(iplon,lay) -1))*nspb(29) + 1
         tauray = colmol(iplon,lay)  * rayl

         do ig = 1, ng29
            taug(iplon,lay,ngs28+ig)  = colco2(iplon,lay)  * &
                (fac00(iplon,lay)  * absb(ind0,ig) + &
                 fac10(iplon,lay)  * absb(ind0+1,ig) + &
                 fac01(iplon,lay)  * absb(ind1,ig) + &
                 fac11(iplon,lay)  * absb(ind1+1,ig)) &  
                 + colh2o(iplon,lay)  * absh2o(ig) 
!            ssa(lay,ngs28+ig) = tauray/taug(lay,ngs28+ig)
         
            taur(iplon,lay,ngs28+ig)  = tauray
         enddo
         end if
         enddo
      enddo

!$acc end kernels
      end subroutine taumol29

     

      end module rrtmg_sw_taumol

