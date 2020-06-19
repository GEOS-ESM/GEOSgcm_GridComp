!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
#include "_gpudef.inc"
#include "_memdef.inc"
      module gpu_rrtmg_lw_taumol

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! (dmb 2012) This is the GPU version of the taumol subroutines.  At first I was going to 
! try and combine the taumol routines into a single subroutine, but it turns out that 
! all 16 can remain and run efficiently on the GPU.  


! ------- Modules -------

       !use parkind, only : im => kind , rb => kind  
      use parrrtm, only : mg, nbndlw, maxxsec, ngptlw
      use rrlw_con, only: oneminus
      use rrlw_wvn, only: nspa, nspb
      use rrlw_vsn, only: hvrtau, hnamtau
      use rrlw_wvn, only: ngb
      use rrlw_ref
      use memory
 

#ifdef _CUDA
      use cudafor
#endif

      implicit none

      ! (dmb 2012) There are a lot of GPU module level variables in this module
      ! The parameter list for the taumol subroutines have been reduced for 
      ! efficiency and readability.
      real  _gpudev, allocatable :: pavel(:,:)
      real  _gpudev, allocatable :: wx1(:,:)
      real  _gpudev, allocatable :: wx2(:,:)
      real  _gpudev, allocatable :: wx3(:,:)
      real  _gpudev, allocatable :: wx4(:,:)
      real  _gpudev, allocatable :: coldry(:,:)
      integer  _gpudev, allocatable :: laytrop(:)
      integer  _gpudev, allocatable :: jp(:,:)
      integer  _gpudev, allocatable :: jt(:,:)
      integer  _gpudev, allocatable :: jt1(:,:)
      real  _gpudev, allocatable :: colh2o(:,:)
      real  _gpudev, allocatable :: colco2(:,:)
      real  _gpudev, allocatable :: colo3(:,:)
      real  _gpudev, allocatable :: coln2o(:,:)
      real  _gpudev, allocatable :: colco(:,:)
      real  _gpudev, allocatable :: colch4(:,:)
      real  _gpudev, allocatable :: colo2(:,:)
      real  _gpudev, allocatable :: colbrd(:,:)
      integer  _gpudev, allocatable :: indself(:,:)
      integer  _gpudev, allocatable :: indfor(:,:)
      real  _gpudev, allocatable :: selffac(:,:)
      real  _gpudev, allocatable :: selffrac(:,:)
      real  _gpudev, allocatable :: forfac(:,:)
      real  _gpudev, allocatable :: forfrac(:,:)
      integer  _gpudev, allocatable :: indminor(:,:)
      real  _gpudev, allocatable :: minorfrac(:,:)
      real  _gpudev, allocatable :: scaleminor(:,:)
      real  _gpudev, allocatable :: scaleminorn2(:,:)
      real  _gpudev, allocatable :: fac00(:,:), fac01(:,:), fac10(:,:), fac11(:,:)
      real  _gpudev, allocatable :: rat_h2oco2(:,:),rat_h2oco2_1(:,:), &
                                            rat_h2oo3(:,:),rat_h2oo3_1(:,:), & !    Dimensions: (nlayers)
                                            rat_h2on2o(:,:),rat_h2on2o_1(:,:), &
                                            rat_h2och4(:,:),rat_h2och4_1(:,:), &
                                            rat_n2oco2(:,:),rat_n2oco2_1(:,:), &
                                            rat_o3co2(:,:),rat_o3co2_1(:,:)

      real  _gpudev, allocatable :: tauaa(:,:,:)
     
      integer  _gpudev, allocatable :: nspad(:)
      integer  _gpudev, allocatable :: nspbd(:)
      real  _gpucon :: oneminusd 


      contains

      !----------------------------------------------------------------------------
      _gpuker subroutine taugb1g( ncol, nlayers, taug, fracsd )

        
!----------------------------------------------------------------------------

! ------- Modifications -------
!  Written by Eli J. Mlawer, Atmospheric & Environmental Research.
!  Revised by Michael J. Iacono, Atmospheric & Environmental Research.
!
!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                          (high key - h2o; high minor - n2)
!
!     note: previous versions of rrtm band 1: 
!           10-250 cm-1 (low - h2o; high - h2o)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng1
      use rrlw_kg01

! ------- Declarations -------

! Local 
      integer  :: lay, ind0, ind1, inds, indf, indm, ig
      real  :: pp, corradj, scalen2, tauself, taufor, taun2
      integer , value, intent(in) :: ncol, nlayers
      real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     

      integer  :: iplon

#ifdef _CUDA
       iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif


! Minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

! Lower atmosphere loop


   if (lay <= laytrop(iplon)) then

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(1) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(1) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)
         pp = pavel(iplon, lay)
         corradj =  1.
         if (pp .lt. 250. ) then
            corradj = 1.  - 0.15  * (250. -pp) / 154.4 
         endif

         scalen2 = colbrd(iplon,lay) * scaleminorn2(iplon,lay)
         do ig = 1, ng1
            tauself = selffac(iplon,lay) * (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor =  forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) -  forrefd(indf,ig))) 
            taun2 = scalen2*(ka_mn2d(indm,ig) + & 
                 minorfrac(iplon,lay) * (ka_mn2d(indm+1,ig) - ka_mn2d(indm,ig)))
            taug(iplon,lay,ig) = corradj * (colh2o(iplon,lay) * &
                (fac00(iplon,lay) * absad(ind0,ig) + &
                 fac10(iplon,lay) * absad(ind0+1,ig) + &
                 fac01(iplon,lay) * absad(ind1,ig) + &
                 fac11(iplon,lay) * absad(ind1+1,ig)) & 
                 + tauself + taufor + taun2)
             fracsd(iplon,lay,ig) = fracrefad(ig)
            
         enddo
   else

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(1) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(1) + 1
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)
         pp = pavel(iplon, lay)
         corradj =  1.  - 0.15  * (pp / 95.6 )

         scalen2 = colbrd(iplon,lay) * scaleminorn2(iplon,lay)
         do ig = 1, ng1
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + &
                 forfrac(iplon,lay) * (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            taun2 = scalen2*(kb_mn2d(indm,ig) + & 
                 minorfrac(iplon,lay) * (kb_mn2d(indm+1,ig) - kb_mn2d(indm,ig)))
            taug(iplon,lay,ig) = corradj * (colh2o(iplon,lay) * &
                (fac00(iplon,lay) * absbd(ind0,ig) + &
                 fac10(iplon,lay) * absbd(ind0+1,ig) + &
                 fac01(iplon,lay) * absbd(ind1,ig) + &
                 fac11(iplon,lay) * absbd(ind1+1,ig)) &  
                 + taufor + taun2)
            fracsd(iplon,lay,ig) = fracrefbd(ig)
         enddo
     endif
#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb1g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb2g( ncol, nlayers , taug, fracsd)
!----------------------------------------------------------------------------
!
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!
!     note: previous version of rrtm band 2: 
!           250 - 500 cm-1 (low - h2o; high - h2o)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng2, ngs1
      use parrrtm, only : ngs1
      use rrlw_kg02

! ------- Declarations -------
     real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
! Local 
      integer  :: lay, ind0, ind1, inds, indf, ig
      real  :: pp, corradj, tauself, taufor
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
      iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
     if (lay <= laytrop(iplon)) then

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(2) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(2) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         pp = pavel(iplon, lay)
         corradj = 1.  - .05  * (pp - 100. ) / 900. 
         do ig = 1, ng2
            tauself = selffac(iplon,lay) * (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor =  forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            taug(iplon,lay,ngs1+ig) = corradj * (colh2o(iplon,lay) * &
                (fac00(iplon,lay) * absad(ind0,ig) + &
                 fac10(iplon,lay) * absad(ind0+1,ig) + &
                 fac01(iplon,lay) * absad(ind1,ig) + &
                 fac11(iplon,lay) * absad(ind1+1,ig)) &
                 + tauself + taufor)
            fracsd(iplon,lay,ngs1+ig) = fracrefad(ig)
         enddo
 else

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(2) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(2) + 1
         indf = indfor(iplon,lay)
         do ig = 1, ng2
            taufor =  forfac(iplon,lay) * (forrefd(indf,ig) + &
                 forfrac(iplon,lay) * (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            taug(iplon,lay,ngs1+ig) = colh2o(iplon,lay) * &
                (fac00(iplon,lay) * absbd(ind0,ig) + &
                 fac10(iplon,lay) * absbd(ind0+1,ig) + &
                 fac01(iplon,lay) * absbd(ind1,ig) + &
                 fac11(iplon,lay) * absbd(ind1+1,ig)) &
                 + taufor
            fracsd(iplon,lay,ngs1+ig) = fracrefbd(ig)
         enddo
   endif
      
#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb2g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb3g( ncol, nlayers, taug, fracsd )
!----------------------------------------------------------------------------
!
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng3, ngs2
      use parrrtm, only : ngs2
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg03

! ------- Declarations -------

! Local 
 real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
      integer  :: lay, ind0, ind1, inds, indf, indm, ig
      integer  :: js, js1, jmn2o, jpl
      real  :: speccomb, specparm, specmult, fs
      real  :: speccomb1, specparm1, specmult1, fs1
      real  :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, &
                       fmn2o, fmn2omf, chi_n2o, ratn2o, adjfac, adjcoln2o
      real  :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real  :: p, p4, fk0, fk1, fk2
      real  :: fac000, fac100, fac200, fac010, fac110, fac210
      real  :: fac001, fac101, fac201, fac011, fac111, fac211
      real  :: tauself, taufor, n2om1, n2om2, absn2o
      real  :: refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b
      real  :: tau_major, tau_major1
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
      iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Minor gas mapping levels:
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

!  P = 212.725 mb
      refrat_planck_a = chi_mlsd(1,9)/chi_mlsd(2,9)

!  P = 95.58 mb
      refrat_planck_b = chi_mlsd(1,13)/chi_mlsd(2,13)

!  P = 706.270mb
      refrat_m_a = chi_mlsd(1,3)/chi_mlsd(2,3)

!  P = 95.58 mb 
      refrat_m_b = chi_mlsd(1,13)/chi_mlsd(2,13)

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature, and appropriate species.  Below laytrop, the water vapor 
! self-continuum and foreign continuum is interpolated (in temperature) 
! separately.

! Lower atmosphere loop
      if (lay <= laytrop(iplon)) then

         speccomb = colh2o(iplon,lay) + rat_h2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )        

         speccomb1 = colh2o(iplon,lay) + rat_h2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mn2o = colh2o(iplon,lay) + refrat_m_a*colco2(iplon,lay)
         specparm_mn2o = colh2o(iplon,lay)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminusd) specparm_mn2o = oneminusd
         specmult_mn2o = 8. *specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0 )
         fmn2omf = minorfrac(iplon,lay)*fmn2o
!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(iplon,lay)/coldry(iplon,lay)
         ratn2o = 1.e20 *chi_n2o/chi_mlsd(4,jp(iplon,lay)+1)
         if (ratn2o .gt. 1.5 ) then
            adjfac = 0.5 +(ratn2o-0.5 )**0.65 
            adjcoln2o = adjfac*chi_mlsd(4,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20 
         else
            adjcoln2o = coln2o(iplon,lay)
         endif

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colco2(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(3) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(3) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1.  - fs) * fac00(iplon,lay)
            fac010 = (1.  - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif
         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1.  - fs1) * fac01(iplon,lay)
            fac011 = (1.  - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng3
            tauself = selffac(iplon,lay)* (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            n2om1 = ka_mn2od(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2od(jmn2o+1,indm,ig) - ka_mn2od(jmn2o,indm,ig))
            n2om2 = ka_mn2od(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2od(jmn2o+1,indm+1,ig) - ka_mn2od(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(iplon,lay) * (n2om2 - n2om1)

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac200 * absad(ind0+2,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig) + &
                    fac210 * absad(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                    (fac200 * absad(ind0-1,ig) + &
                    fac100 * absad(ind0,ig) + &
                    fac000 * absad(ind0+1,ig) + &
                    fac210 * absad(ind0+8,ig) + &
                    fac110 * absad(ind0+9,ig) + &
                    fac010 * absad(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac201 * absad(ind1+2,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig) + &
                    fac211 * absad(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                    (fac201 * absad(ind1-1,ig) + &
                    fac101 * absad(ind1,ig) + &
                    fac001 * absad(ind1+1,ig) + &
                    fac211 * absad(ind1+8,ig) + &
                    fac111 * absad(ind1+9,ig) + &
                    fac011 * absad(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) +  &
                    fac101 * absad(ind1+1,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig))
            endif

            taug(iplon,lay,ngs2+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracsd(iplon,lay,ngs2+ig) = fracrefad(ig,jpl) + fpl * &
                 (fracrefad(ig,jpl+1)-fracrefad(ig,jpl))
         enddo
    

! Upper atmosphere loop
     else

         speccomb = colh2o(iplon,lay) + rat_h2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 4. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(iplon,lay) + rat_h2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 4. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         fac000 = (1.  - fs) * fac00(iplon,lay)
         fac010 = (1.  - fs) * fac10(iplon,lay)
         fac100 = fs * fac00(iplon,lay)
         fac110 = fs * fac10(iplon,lay)
         fac001 = (1.  - fs1) * fac01(iplon,lay)
         fac011 = (1.  - fs1) * fac11(iplon,lay)
         fac101 = fs1 * fac01(iplon,lay)
         fac111 = fs1 * fac11(iplon,lay)

         speccomb_mn2o = colh2o(iplon,lay) + refrat_m_b*colco2(iplon,lay)
         specparm_mn2o = colh2o(iplon,lay)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminusd) specparm_mn2o = oneminusd
         specmult_mn2o = 4. *specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0 )
         fmn2omf = minorfrac(iplon,lay)*fmn2o
!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(iplon,lay)/coldry(iplon,lay)
         ratn2o = 1.e20*chi_n2o/chi_mlsd(4,jp(iplon,lay)+1)
         if (ratn2o .gt. 1.5 ) then
            adjfac = 0.5 +(ratn2o-0.5 )**0.65 
            adjcoln2o = adjfac*chi_mlsd(4,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20 
         else
            adjcoln2o = coln2o(iplon,lay)
         endif

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_b*colco2(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 4. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(3) + js
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(3) + js1
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         do ig = 1, ng3
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + &
                 forfrac(iplon,lay) * (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            n2om1 = kb_mn2od(jmn2o,indm,ig) + fmn2o * &
                 (kb_mn2od(jmn2o+1,indm,ig)-kb_mn2od(jmn2o,indm,ig))
            n2om2 = kb_mn2od(jmn2o,indm+1,ig) + fmn2o * &
                 (kb_mn2od(jmn2o+1,indm+1,ig)-kb_mn2od(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(iplon,lay) * (n2om2 - n2om1)
            taug(iplon,lay,ngs2+ig) = speccomb * &
                (fac000 * absbd(ind0,ig) + &
                fac100 * absbd(ind0+1,ig) + &
                fac010 * absbd(ind0+5,ig) + &
                fac110 * absbd(ind0+6,ig)) &
                + speccomb1 * &
                (fac001 * absbd(ind1,ig) +  &
                fac101 * absbd(ind1+1,ig) + &
                fac011 * absbd(ind1+5,ig) + &
                fac111 * absbd(ind1+6,ig))  &
                + taufor &
                + adjcoln2o*absn2o
            fracsd(iplon,lay,ngs2+ig) = fracrefbd(ig,jpl) + fpl * &
                (fracrefbd(ig,jpl+1)-fracrefbd(ig,jpl))
         enddo
     endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb3g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb4g( ncol, nlayers, taug, fracsd )
!----------------------------------------------------------------------------
!
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng4, ngs3
      use parrrtm, only : ngs3
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg04

! ------- Declarations -------

! Local 
       real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
      
      integer  :: lay, ind0, ind1, inds, indf, ig
      integer  :: js, js1, jpl
      real  :: speccomb, specparm, specmult, fs
      real  :: speccomb1, specparm1, specmult1, fs1
      real  :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real  :: p, p4, fk0, fk1, fk2
      real  :: fac000, fac100, fac200, fac010, fac110, fac210
      real  :: fac001, fac101, fac201, fac011, fac111, fac211
      real  :: tauself, taufor
      real  :: refrat_planck_a, refrat_planck_b
      real  :: tau_major, tau_major1
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
        iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif 
! P =   142.5940 mb
      refrat_planck_a = chi_mlsd(1,11)/chi_mlsd(2,11)

! P = 95.58350 mb
      refrat_planck_b = chi_mlsd(3,13)/chi_mlsd(2,13)

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature, and appropriate species.  Below laytrop, the water 
! vapor self-continuum and foreign continuum is interpolated (in temperature) 
! separately.

! Lower atmosphere loop
     if (lay <= laytrop(iplon)) then

         speccomb = colh2o(iplon,lay) + rat_h2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(iplon,lay) + rat_h2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colco2(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(4) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(4) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1.  - fs) * fac00(iplon,lay)
            fac010 = (1.  - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1.  - fs1) * fac01(iplon,lay)
            fac011 = (1.  - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng4
            tauself = selffac(iplon,lay)* (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor =  forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac200 * absad(ind0+2,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig) + &
                    fac210 * absad(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                    (fac200 * absad(ind0-1,ig) + &
                    fac100 * absad(ind0,ig) + &
                    fac000 * absad(ind0+1,ig) + &
                    fac210 * absad(ind0+8,ig) + &
                    fac110 * absad(ind0+9,ig) + &
                    fac010 * absad(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) +  &
                    fac101 * absad(ind1+1,ig) + &
                    fac201 * absad(ind1+2,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig) + &
                    fac211 * absad(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                    (fac201 * absad(ind1-1,ig) + &
                    fac101 * absad(ind1,ig) + &
                    fac001 * absad(ind1+1,ig) + &
                    fac211 * absad(ind1+8,ig) + &
                    fac111 * absad(ind1+9,ig) + &
                    fac011 * absad(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig))
            endif

            taug(iplon,lay,ngs3+ig) = tau_major + tau_major1 &
                 + tauself + taufor
            fracsd(iplon,lay,ngs3+ig) = fracrefad(ig,jpl) + fpl * &
                 (fracrefad(ig,jpl+1)-fracrefad(ig,jpl))
         enddo
    

! Upper atmosphere loop
     else

         speccomb = colo3(iplon,lay) + rat_o3co2(iplon,lay)*colco2(iplon,lay)
         specparm = colo3(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 4. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colo3(iplon,lay) + rat_o3co2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colo3(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 4. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         fac000 = (1.  - fs) * fac00(iplon,lay)
         fac010 = (1.  - fs) * fac10(iplon,lay)
         fac100 = fs * fac00(iplon,lay)
         fac110 = fs * fac10(iplon,lay)
         fac001 = (1.  - fs1) * fac01(iplon,lay)
         fac011 = (1.  - fs1) * fac11(iplon,lay)
         fac101 = fs1 * fac01(iplon,lay)
         fac111 = fs1 * fac11(iplon,lay)

         speccomb_planck = colo3(iplon,lay)+refrat_planck_b*colco2(iplon,lay)
         specparm_planck = colo3(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 4. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(4) + js
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(4) + js1

         do ig = 1, ng4
            taug(iplon,lay,ngs3+ig) =  speccomb * &
                (fac000 * absbd(ind0,ig) + &
                fac100 * absbd(ind0+1,ig) + &
                fac010 * absbd(ind0+5,ig) + &
                fac110 * absbd(ind0+6,ig)) &
                + speccomb1 * &
                (fac001 * absbd(ind1,ig) +  &
                fac101 * absbd(ind1+1,ig) + &
                fac011 * absbd(ind1+5,ig) + &
                fac111 * absbd(ind1+6,ig))
            fracsd(iplon,lay,ngs3+ig) = fracrefbd(ig,jpl) + fpl * &
                (fracrefbd(ig,jpl+1)-fracrefbd(ig,jpl))
         enddo

! Empirical modification to code to improve stratospheric cooling rates
! for co2.  Revised to apply weighting for g-point reduction in this band.

         taug(iplon,lay,ngs3+8)=taug(iplon,lay,ngs3+8)*0.92
         taug(iplon,lay,ngs3+9)=taug(iplon,lay,ngs3+9)*0.88
         taug(iplon,lay,ngs3+10)=taug(iplon,lay,ngs3+10)*1.07
         taug(iplon,lay,ngs3+11)=taug(iplon,lay,ngs3+11)*1.1
         taug(iplon,lay,ngs3+12)=taug(iplon,lay,ngs3+12)*0.99
         taug(iplon,lay,ngs3+13)=taug(iplon,lay,ngs3+13)*0.88
         taug(iplon,lay,ngs3+14)=taug(iplon,lay,ngs3+14)*0.943

     endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb4g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb5g( ncol, nlayers , taug, fracsd)
!----------------------------------------------------------------------------
!
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                           (high key - o3,co2)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng5, ngs4
      use parrrtm, only : ngs4
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg05

! ------- Declarations -------

! Local 
      real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
      integer  :: lay, ind0, ind1, inds, indf, indm, ig
      integer  :: js, js1, jmo3, jpl
      real  :: speccomb, specparm, specmult, fs
      real  :: speccomb1, specparm1, specmult1, fs1
      real  :: speccomb_mo3, specparm_mo3, specmult_mo3, fmo3
      real  :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real  :: p, p4, fk0, fk1, fk2
      real  :: fac000, fac100, fac200, fac010, fac110, fac210
      real  :: fac001, fac101, fac201, fac011, fac111, fac211
      real  :: tauself, taufor, o3m1, o3m2, abso3
      real  :: refrat_planck_a, refrat_planck_b, refrat_m_a
      real  :: tau_major, tau_major1
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
      iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Minor gas mapping level :
!     lower - o3, p = 317.34 mbar, t = 240.77 k
!     lower - ccl4

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

! P = 473.420 mb
      refrat_planck_a = chi_mlsd(1,5)/chi_mlsd(2,5)

! P = 0.2369 mb
      refrat_planck_b = chi_mlsd(3,43)/chi_mlsd(2,43)

! P = 317.3480
      refrat_m_a = chi_mlsd(1,7)/chi_mlsd(2,7)

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature, and appropriate species.  Below laytrop, the 
! water vapor self-continuum and foreign continuum is 
! interpolated (in temperature) separately.

! Lower atmosphere loop
      !do lay = 1, laytrop(iplon)
      if (lay <= laytrop(iplon)) then
         speccomb = colh2o(iplon,lay) + rat_h2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(iplon,lay) + rat_h2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mo3 = colh2o(iplon,lay) + refrat_m_a*colco2(iplon,lay)
         specparm_mo3 = colh2o(iplon,lay)/speccomb_mo3
         if (specparm_mo3 .ge. oneminusd) specparm_mo3 = oneminusd
         specmult_mo3 = 8. *specparm_mo3
         jmo3 = 1 + int(specmult_mo3)
         fmo3 = mod(specmult_mo3,1.0 )

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colco2(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(5) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(5) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1.  - fs) * fac00(iplon,lay)
            fac010 = (1.  - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1.  - fs1) * fac01(iplon,lay)
            fac011 = (1.  - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng5
            tauself = selffac(iplon,lay) * (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor =  forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            o3m1 = ka_mo3d(jmo3,indm,ig) + fmo3 * &
                 (ka_mo3d(jmo3+1,indm,ig)-ka_mo3d(jmo3,indm,ig))
            o3m2 = ka_mo3d(jmo3,indm+1,ig) + fmo3 * &
                 (ka_mo3d(jmo3+1,indm+1,ig)-ka_mo3d(jmo3,indm+1,ig))
            abso3 = o3m1 + minorfrac(iplon,lay)*(o3m2-o3m1)

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac200 * absad(ind0+2,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig) + &
                    fac210 * absad(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                    (fac200 * absad(ind0-1,ig) + &
                    fac100 * absad(ind0,ig) + &
                    fac000 * absad(ind0+1,ig) + &
                    fac210 * absad(ind0+8,ig) + &
                    fac110 * absad(ind0+9,ig) + &
                    fac010 * absad(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac201 * absad(ind1+2,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig) + &
                    fac211 * absad(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * & 
                    (fac201 * absad(ind1-1,ig) + &
                    fac101 * absad(ind1,ig) + &
                    fac001 * absad(ind1+1,ig) + &
                    fac211 * absad(ind1+8,ig) + &
                    fac111 * absad(ind1+9,ig) + &
                    fac011 * absad(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig))
            endif

            taug(iplon,lay,ngs4+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + abso3*colo3(iplon,lay) &
                 + wx1(iplon,lay) * coldry(iplon,lay) * 1.e-20  * ccl4d(ig)
            fracsd(iplon,lay,ngs4+ig) = fracrefad(ig,jpl) + fpl * &
                 (fracrefad(ig,jpl+1)-fracrefad(ig,jpl))
         enddo
      else

! Upper atmosphere loop
      !do lay = laytrop(iplon)+1, nlayers

         speccomb = colo3(iplon,lay) + rat_o3co2(iplon,lay)*colco2(iplon,lay)
         specparm = colo3(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 4. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colo3(iplon,lay) + rat_o3co2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colo3(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 4. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         fac000 = (1.  - fs) * fac00(iplon,lay)
         fac010 = (1.  - fs) * fac10(iplon,lay)
         fac100 = fs * fac00(iplon,lay)
         fac110 = fs * fac10(iplon,lay)
         fac001 = (1.  - fs1) * fac01(iplon,lay)
         fac011 = (1.  - fs1) * fac11(iplon,lay)
         fac101 = fs1 * fac01(iplon,lay)
         fac111 = fs1 * fac11(iplon,lay)

         speccomb_planck = colo3(iplon,lay)+refrat_planck_b*colco2(iplon,lay)
         specparm_planck = colo3(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 4. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(5) + js
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(5) + js1
         
         do ig = 1, ng5
            taug(iplon,lay,ngs4+ig) = speccomb * &
                (fac000 * absbd(ind0,ig) + &
                fac100 * absbd(ind0+1,ig) + &
                fac010 * absbd(ind0+5,ig) + &
                fac110 * absbd(ind0+6,ig)) &
                + speccomb1 * &
                (fac001 * absbd(ind1,ig) + &
                fac101 * absbd(ind1+1,ig) + &
                fac011 * absbd(ind1+5,ig) + &
                fac111 * absbd(ind1+6,ig))  &
                + wx1(iplon, lay) * coldry(iplon,lay) * 1.e-20  * ccl4d(ig)
            fracsd(iplon,lay,ngs4+ig) = fracrefbd(ig,jpl) + fpl * &
                (fracrefbd(ig,jpl+1)-fracrefbd(ig,jpl))
         enddo
      endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb5g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb6g( ncol, nlayers, taug, fracsd )
!----------------------------------------------------------------------------
!
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng6, ngs5
      use parrrtm, only : ngs5
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg06

! ------- Declarations -------

! Local 
      integer  :: lay, ind0, ind1, inds, indf, indm, ig
      real  :: chi_co2, ratco2, adjfac, adjcolco2
      real  :: tauself, taufor, absco2
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
      real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
#ifdef _CUDA
        iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12

! Compute the optical depth by interpolating in ln(pressure) and
! temperature. The water vapor self-continuum and foreign continuum
! is interpolated (in temperature) separately.  

! Lower atmosphere loop
     if (lay <= laytrop(iplon)) then

! In atmospheres where the amount of CO2 is too great to be considered
! a minor species, adjust the column amount of CO2 by an empirical factor 
! to obtain the proper contribution.
         chi_co2 = colco2(iplon,lay)/(coldry(iplon,lay))
         ratco2 = 1.e20 *chi_co2/chi_mlsd(2,jp(iplon,lay)+1)
         if (ratco2 .gt. 3.0 ) then
            adjfac = 2.0 +(ratco2-2.0 )**0.77 
            adjcolco2 = adjfac*chi_mlsd(2,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20 
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(6) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(6) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         do ig = 1, ng6
            tauself = selffac(iplon,lay) * (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor =  forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig)))
            absco2 =  (ka_mco2d(indm,ig) + minorfrac(iplon,lay) * &
                 (ka_mco2d(indm+1,ig) - ka_mco2d(indm,ig)))
            taug(iplon,lay,ngs5+ig) = colh2o(iplon,lay) * &
                (fac00(iplon,lay) * absad(ind0,ig) + &
                 fac10(iplon,lay) * absad(ind0+1,ig) + &
                 fac01(iplon,lay) * absad(ind1,ig) +  &
                 fac11(iplon,lay) * absad(ind1+1,ig))  &
                 + tauself + taufor &
                 + adjcolco2 * absco2 &
                 + wx2(iplon, lay) * coldry(iplon,lay) * 1.e-20  * cfc11adjd(ig) &
                 + wx3(iplon, lay) * coldry(iplon,lay) * 1.e-20  * cfc12d(ig)
            fracsd(iplon,lay,ngs5+ig) = fracrefad(ig)
         enddo
   else

         do ig = 1, ng6
            taug(iplon,lay,ngs5+ig) = 0.0  &
                 + wx2(iplon, lay) * coldry(iplon,lay) * 1.e-20   * cfc11adjd(ig) &
                 + wx3(iplon, lay) * coldry(iplon,lay) * 1.e-20  * cfc12d(ig)
            fracsd(iplon,lay,ngs5+ig) = fracrefad(ig)
         enddo
     endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb6g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb7g( ncol, nlayers , taug, fracsd)
!----------------------------------------------------------------------------
!
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!                            (high key - o3; high minor - co2)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng7, ngs6
      use parrrtm, only : ngs6
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg07

! ------- Declarations -------

! Local 
      real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
      integer  :: lay, ind0, ind1, inds, indf, indm, ig
      integer  :: js, js1, jmco2, jpl
      real  :: speccomb, specparm, specmult, fs
      real  :: speccomb1, specparm1, specmult1, fs1
      real  :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real  :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real  :: p, p4, fk0, fk1, fk2
      real  :: fac000, fac100, fac200, fac010, fac110, fac210
      real  :: fac001, fac101, fac201, fac011, fac111, fac211
      real  :: tauself, taufor, co2m1, co2m2, absco2
      real  :: chi_co2, ratco2, adjfac, adjcolco2
      real  :: refrat_planck_a, refrat_m_a
      real  :: tau_major, tau_major1
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
      iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Minor gas mapping level :
!     lower - co2, p = 706.2620 mbar, t= 278.94 k
!     upper - co2, p = 12.9350 mbar, t = 234.01 k

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower atmosphere.

! P = 706.2620 mb
      refrat_planck_a = chi_mlsd(1,3)/chi_mlsd(3,3)

! P = 706.2720 mb
      refrat_m_a = chi_mlsd(1,3)/chi_mlsd(3,3)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately. 

! Lower atmosphere loop
      if (lay <= laytrop(iplon)) then

         speccomb = colh2o(iplon,lay) + rat_h2oo3(iplon,lay)*colo3(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(iplon,lay) + rat_h2oo3_1(iplon,lay)*colo3(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mco2 = colh2o(iplon,lay) + refrat_m_a*colo3(iplon,lay)
         specparm_mco2 = colh2o(iplon,lay)/speccomb_mco2
         if (specparm_mco2 .ge. oneminusd) specparm_mco2 = oneminusd
         specmult_mco2 = 8. *specparm_mco2

         jmco2 = 1 + int(specmult_mco2)
         fmco2 = mod(specmult_mco2,1.0 )

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(iplon,lay)/(coldry(iplon,lay))
         ratco2 = 1.e20*chi_co2/chi_mlsd(2,jp(iplon,lay)+1)
         if (ratco2 .gt. 3.0 ) then
            adjfac = 3.0 +(ratco2-3.0 )**0.79 
            adjcolco2 = adjfac*chi_mlsd(2,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20 
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colo3(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(7) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(7) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1.  - fs) * fac00(iplon,lay)
            fac010 = (1.  - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif
         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1.  - fs1) * fac01(iplon,lay)
            fac011 = (1.  - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng7
            tauself = selffac(iplon,lay)* (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            co2m1 = ka_mco2d(jmco2,indm,ig) + fmco2 * &
                 (ka_mco2d(jmco2+1,indm,ig) - ka_mco2d(jmco2,indm,ig))
            co2m2 = ka_mco2d(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco2d(jmco2+1,indm+1,ig) - ka_mco2d(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(iplon,lay) * (co2m2 - co2m1)

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac200 * absad(ind0+2,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig) + &
                    fac210 * absad(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                    (fac200 * absad(ind0-1,ig) + &
                    fac100 * absad(ind0,ig) + &
                    fac000 * absad(ind0+1,ig) + &
                    fac210 * absad(ind0+8,ig) + &
                    fac110 * absad(ind0+9,ig) + &
                    fac010 * absad(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac201 * absad(ind1+2,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig) + &
                    fac211 * absad(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                    (fac201 * absad(ind1-1,ig) + &
                    fac101 * absad(ind1,ig) + &
                    fac001 * absad(ind1+1,ig) + &
                    fac211 * absad(ind1+8,ig) + &
                    fac111 * absad(ind1+9,ig) + &
                    fac011 * absad(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) +  &
                    fac101 * absad(ind1+1,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig))
            endif

            taug(iplon,lay,ngs6+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcolco2*absco2
            fracsd(iplon,lay,ngs6+ig) = fracrefad(ig,jpl) + fpl * &
                 (fracrefad(ig,jpl+1)-fracrefad(ig,jpl))
         enddo
    else
!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(iplon,lay)/(coldry(iplon,lay))
         ratco2 = 1.e20*chi_co2/chi_mlsd(2,jp(iplon,lay)+1)
         if (ratco2 .gt. 3.0 ) then
            adjfac = 2.0 +(ratco2-2.0 )**0.79 
            adjcolco2 = adjfac*chi_mlsd(2,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20 
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(7) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(7) + 1
         indm = indminor(iplon,lay)

         do ig = 1, ng7
            absco2 = kb_mco2d(indm,ig) + minorfrac(iplon,lay) * &
                 (kb_mco2d(indm+1,ig) - kb_mco2d(indm,ig))
            taug(iplon,lay,ngs6+ig) = colo3(iplon,lay) * &
                 (fac00(iplon,lay) * absbd(ind0,ig) + &
                 fac10(iplon,lay) * absbd(ind0+1,ig) + &
                 fac01(iplon,lay) * absbd(ind1,ig) + &
                 fac11(iplon,lay) * absbd(ind1+1,ig)) &
                 + adjcolco2 * absco2
            fracsd(iplon,lay,ngs6+ig) = fracrefbd(ig)
         enddo

! Empirical modification to code to improve stratospheric cooling rates
! for o3.  Revised to apply weighting for g-point reduction in this band.

         taug(iplon,lay,ngs6+6)=taug(iplon,lay,ngs6+6)*0.92 
         taug(iplon,lay,ngs6+7)=taug(iplon,lay,ngs6+7)*0.88 
         taug(iplon,lay,ngs6+8)=taug(iplon,lay,ngs6+8)*1.07 
         taug(iplon,lay,ngs6+9)=taug(iplon,lay,ngs6+9)*1.1 
         taug(iplon,lay,ngs6+10)=taug(iplon,lay,ngs6+10)*0.99 
         taug(iplon,lay,ngs6+11)=taug(iplon,lay,ngs6+11)*0.855 

      endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb7g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb8g( ncol, nlayers, taug, fracsd )
!----------------------------------------------------------------------------
!
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                             (high key - o3; high minor - co2, n2o)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng8, ngs7
      use parrrtm, only : ngs7
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg08

! ------- Declarations -------

! Local  
        real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
      integer  :: lay, ind0, ind1, inds, indf, indm, ig
      real  :: tauself, taufor, absco2, abso3, absn2o
      real  :: chi_co2, ratco2, adjfac, adjcolco2
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
      iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature, and appropriate species.  Below laytrop, the water vapor 
! self-continuum and foreign continuum is interpolated (in temperature) 
! separately.

! Lower atmosphere loop
     if (lay <= laytrop(iplon)) then

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(iplon,lay)/(coldry(iplon,lay))
         ratco2 = 1.e20 *chi_co2/chi_mlsd(2,jp(iplon,lay)+1)
         if (ratco2 .gt. 3.0 ) then
            adjfac = 2.0 +(ratco2-2.0 )**0.65 
            adjcolco2 = adjfac*chi_mlsd(2,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20 
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(8) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(8) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         do ig = 1, ng8
            tauself = selffac(iplon,lay) * (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig)))
            absco2 =  (ka_mco2d(indm,ig) + minorfrac(iplon,lay) * &
                 (ka_mco2d(indm+1,ig) - ka_mco2d(indm,ig)))
            abso3 =  (ka_mo3d(indm,ig) + minorfrac(iplon,lay) * &
                 (ka_mo3d(indm+1,ig) - ka_mo3d(indm,ig)))
            absn2o =  (ka_mn2od(indm,ig) + minorfrac(iplon,lay) * &
                 (ka_mn2od(indm+1,ig) - ka_mn2od(indm,ig)))
            taug(iplon,lay,ngs7+ig) = colh2o(iplon,lay) * &
                 (fac00(iplon,lay) * absad(ind0,ig) + &
                 fac10(iplon,lay) * absad(ind0+1,ig) + &
                 fac01(iplon,lay) * absad(ind1,ig) +  &
                 fac11(iplon,lay) * absad(ind1+1,ig)) &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colo3(iplon,lay) * abso3 &
                 + coln2o(iplon,lay) * absn2o &
                 + wx3(iplon, lay) * coldry(iplon,lay) * 1.e-20  * cfc12d(ig) &
                 + wx4(iplon, lay) * coldry(iplon,lay) * 1.e-20  * cfc22adjd(ig)
            fracsd(iplon,lay,ngs7+ig) = fracrefad(ig)
         enddo
  else
!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(iplon,lay)/coldry(iplon,lay)
         ratco2 = 1.e20 *chi_co2/chi_mlsd(2,jp(iplon,lay)+1)
         if (ratco2 .gt. 3.0 ) then
            adjfac = 2.0 +(ratco2-2.0 )**0.65 
            adjcolco2 = adjfac*chi_mlsd(2,jp(iplon,lay)+1) * coldry(iplon,lay)*1.e-20 
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(8) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(8) + 1
         indm = indminor(iplon,lay)

         do ig = 1, ng8
            absco2 =  (kb_mco2d(indm,ig) + minorfrac(iplon,lay) * &
                 (kb_mco2d(indm+1,ig) - kb_mco2d(indm,ig)))
            absn2o =  (kb_mn2od(indm,ig) + minorfrac(iplon,lay) * &
                 (kb_mn2od(indm+1,ig) - kb_mn2od(indm,ig)))
            taug(iplon,lay,ngs7+ig) = colo3(iplon,lay) * &
                 (fac00(iplon,lay) * absbd(ind0,ig) + &
                 fac10(iplon,lay) * absbd(ind0+1,ig) + &
                 fac01(iplon,lay) * absbd(ind1,ig) + &
                 fac11(iplon,lay) * absbd(ind1+1,ig)) &
                 + adjcolco2*absco2 &
                 + coln2o(iplon,lay)*absn2o & 
                 + wx3(iplon,lay) * coldry(iplon,lay) * 1.e-20  * cfc12d(ig) &
                 + wx4(iplon,lay) * coldry(iplon,lay) * 1.e-20  * cfc22adjd(ig)
            fracsd(iplon,lay,ngs7+ig) = fracrefbd(ig)
         enddo
      endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb8g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb9g( ncol, nlayers, taug, fracsd )
!----------------------------------------------------------------------------
!
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                             (high key - ch4; high minor - n2o)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng9, ngs8
      use parrrtm, only : ngs8
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg09

! ------- Declarations -------
 real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
! Local 
      integer  :: lay, ind0, ind1, inds, indf, indm, ig
      integer  :: js, js1, jmn2o, jpl
      real  :: speccomb, specparm, specmult, fs
      real  :: speccomb1, specparm1, specmult1, fs1
      real  :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o
      real  :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real  :: p, p4, fk0, fk1, fk2
      real  :: fac000, fac100, fac200, fac010, fac110, fac210
      real  :: fac001, fac101, fac201, fac011, fac111, fac211
      real  :: tauself, taufor, n2om1, n2om2, absn2o
      real  :: chi_n2o, ratn2o, adjfac, adjcoln2o
      real  :: refrat_planck_a, refrat_m_a
      real  :: tau_major, tau_major1
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
        iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Minor gas mapping level :
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

! P = 212 mb
      refrat_planck_a = chi_mlsd(1,9)/chi_mlsd(6,9)

! P = 706.272 mb 
      refrat_m_a = chi_mlsd(1,3)/chi_mlsd(6,3)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately.  

! Lower atmosphere loop
      if (lay <= laytrop(iplon)) then

         speccomb = colh2o(iplon,lay) + rat_h2och4(iplon,lay)*colch4(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(iplon,lay) + rat_h2och4_1(iplon,lay)*colch4(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mn2o = colh2o(iplon,lay) + refrat_m_a*colch4(iplon,lay)
         specparm_mn2o = colh2o(iplon,lay)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminusd) specparm_mn2o = oneminusd
         specmult_mn2o = 8. *specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0 )

!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(iplon,lay)/(coldry(iplon,lay))
         ratn2o = 1.e20 *chi_n2o/chi_mlsd(4,jp(iplon,lay)+1)
         if (ratn2o .gt. 1.5 ) then
            adjfac = 0.5 +(ratn2o-0.5 )**0.65 
            adjcoln2o = adjfac*chi_mlsd(4,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20 
         else
            adjcoln2o = coln2o(iplon,lay)
         endif

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colch4(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(9) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(9) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1.  - fs) * fac00(iplon,lay)
            fac010 = (1.  - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1.  - fs1) * fac01(iplon,lay)
            fac011 = (1.  - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng9
            tauself = selffac(iplon,lay)* (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            n2om1 = ka_mn2od(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2od(jmn2o+1,indm,ig) - ka_mn2od(jmn2o,indm,ig))
            n2om2 = ka_mn2od(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2od(jmn2o+1,indm+1,ig) - ka_mn2od(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(iplon,lay) * (n2om2 - n2om1)

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac200 * absad(ind0+2,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig) + &
                    fac210 * absad(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                    (fac200 * absad(ind0-1,ig) + &
                    fac100 * absad(ind0,ig) + &
                    fac000 * absad(ind0+1,ig) + &
                    fac210 * absad(ind0+8,ig) + &
                    fac110 * absad(ind0+9,ig) + &
                    fac010 * absad(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + & 
                    fac101 * absad(ind1+1,ig) + &
                    fac201 * absad(ind1+2,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig) + &
                    fac211 * absad(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                    (fac201 * absad(ind1-1,ig) + &
                    fac101 * absad(ind1,ig) + &
                    fac001 * absad(ind1+1,ig) + &
                    fac211 * absad(ind1+8,ig) + &
                    fac111 * absad(ind1+9,ig) + &
                    fac011 * absad(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig))
            endif

            taug(iplon,lay,ngs8+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracsd(iplon,lay,ngs8+ig) = fracrefad(ig,jpl) + fpl * &
                 (fracrefad(ig,jpl+1)-fracrefad(ig,jpl))
         enddo
 else
!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(iplon,lay)/(coldry(iplon,lay))
         ratn2o = 1.e20 *chi_n2o/chi_mlsd(4,jp(iplon,lay)+1)
         if (ratn2o .gt. 1.5 ) then
            adjfac = 0.5 +(ratn2o-0.5 )**0.65 
            adjcoln2o = adjfac*chi_mlsd(4,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20 
         else
            adjcoln2o = coln2o(iplon,lay)
         endif

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(9) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(9) + 1
         indm = indminor(iplon,lay)

         do ig = 1, ng9
            absn2o = kb_mn2od(indm,ig) + minorfrac(iplon,lay) * &
                (kb_mn2od(indm+1,ig) - kb_mn2od(indm,ig))
            taug(iplon,lay,ngs8+ig) = colch4(iplon,lay) * &
                 (fac00(iplon,lay) * absbd(ind0,ig) + &
                 fac10(iplon,lay) * absbd(ind0+1,ig) + &
                 fac01(iplon,lay) * absbd(ind1,ig) +  &
                 fac11(iplon,lay) * absbd(ind1+1,ig)) &
                 + adjcoln2o*absn2o
            fracsd(iplon,lay,ngs8+ig) = fracrefbd(ig)
         enddo
      endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb9g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb10g( ncol, nlayers, taug, fracsd )
!----------------------------------------------------------------------------
!
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng10, ngs9
      use parrrtm, only : ngs9
      use rrlw_kg10

! ------- Declarations -------
 real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
! Local 
      integer  :: lay, ind0, ind1, inds, indf, ig
      real  :: tauself, taufor
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
     
      iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
     if (lay <= laytrop(iplon)) then
         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(10) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(10) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)

         do ig = 1, ng10
            tauself = selffac(iplon,lay) * (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            taug(iplon,lay,ngs9+ig) = colh2o(iplon,lay) * &
                 (fac00(iplon,lay) * absad(ind0,ig) + &
                 fac10(iplon,lay) * absad(ind0+1,ig) + &
                 fac01(iplon,lay) * absad(ind1,ig) + &
                 fac11(iplon,lay) * absad(ind1+1,ig))  &
                 + tauself + taufor
            fracsd(iplon,lay,ngs9+ig) = fracrefad(ig)
         enddo
   else
   
         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(10) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(10) + 1
         indf = indfor(iplon,lay)

         do ig = 1, ng10
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            taug(iplon,lay,ngs9+ig) = colh2o(iplon,lay) * &
                 (fac00(iplon,lay) * absbd(ind0,ig) + &
                 fac10(iplon,lay) * absbd(ind0+1,ig) + &
                 fac01(iplon,lay) * absbd(ind1,ig) +  &
                 fac11(iplon,lay) * absbd(ind1+1,ig)) &
                 + taufor
            fracsd(iplon,lay,ngs9+ig) = fracrefbd(ig)
         enddo
    end if

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif
      end subroutine taugb10g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb11g( ncol, nlayers, taug, fracsd )
!----------------------------------------------------------------------------
!
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng11, ngs10
      use parrrtm, only : ngs10
      use rrlw_kg11

! ------- Declarations -------
 real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
! Local 
      integer  :: lay, ind0, ind1, inds, indf, indm, ig
      real  :: scaleo2, tauself, taufor, tauo2
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
         iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Minor gas mapping level :
!     lower - o2, p = 706.2720 mbar, t = 278.94 k
!     upper - o2, p = 4.758820 mbarm t = 250.85 k

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

! Lower atmosphere loop
  if (lay <= laytrop(iplon)) then
         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(11) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(11) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)
         scaleo2 = colo2(iplon,lay)*scaleminor(iplon,lay)
         do ig = 1, ng11
            tauself = selffac(iplon,lay) * (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig)))
            tauo2 =  scaleo2 * (ka_mo2d(indm,ig) + minorfrac(iplon,lay) * &
                 (ka_mo2d(indm+1,ig) - ka_mo2d(indm,ig)))
            taug(iplon,lay,ngs10+ig) = colh2o(iplon,lay) * &
                 (fac00(iplon,lay) * absad(ind0,ig) + &
                 fac10(iplon,lay) * absad(ind0+1,ig) + &
                 fac01(iplon,lay) * absad(ind1,ig) + &
                 fac11(iplon,lay) * absad(ind1+1,ig)) &
                 + tauself + taufor &
                 + tauo2
            fracsd(iplon,lay,ngs10+ig) = fracrefad(ig)
         enddo
   else
         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(11) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(11) + 1
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)
         scaleo2 = colo2(iplon,lay)*scaleminor(iplon,lay)
         do ig = 1, ng11
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            tauo2 =  scaleo2 * (kb_mo2d(indm,ig) + minorfrac(iplon,lay) * &
                 (kb_mo2d(indm+1,ig) - kb_mo2d(indm,ig)))
            taug(iplon,lay,ngs10+ig) = colh2o(iplon,lay) * &
                 (fac00(iplon,lay) * absbd(ind0,ig) + &
                 fac10(iplon,lay) * absbd(ind0+1,ig) + &
                 fac01(iplon,lay) * absbd(ind1,ig) + &
                 fac11(iplon,lay) * absbd(ind1+1,ig))  &
                 + taufor &
                 + tauo2
            fracsd(iplon,lay,ngs10+ig) = fracrefbd(ig)
         enddo
      endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb11g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb12g( ncol, nlayers, taug, fracsd )
!----------------------------------------------------------------------------
!
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng12, ngs11
      use parrrtm, only : ngs11
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg12

! ------- Declarations -------
 real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
! Local 
      integer  :: lay, ind0, ind1, inds, indf, ig
      integer  :: js, js1, jpl
      real  :: speccomb, specparm, specmult, fs
      real  :: speccomb1, specparm1, specmult1, fs1
      real  :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real  :: p, p4, fk0, fk1, fk2
      real  :: fac000, fac100, fac200, fac010, fac110, fac210
      real  :: fac001, fac101, fac201, fac011, fac111, fac211
      real  :: tauself, taufor
      real  :: refrat_planck_a
      real  :: tau_major, tau_major1
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
     iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

! P =   174.164 mb 
      refrat_planck_a = chi_mlsd(1,10)/chi_mlsd(2,10)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum adn foreign continuum is interpolated 
! (in temperature) separately.  

! Lower atmosphere loop
    if (lay <= laytrop(iplon)) then

         speccomb = colh2o(iplon,lay) + rat_h2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(iplon,lay) + rat_h2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colco2(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(12) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(12) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1.  - fs) * fac00(iplon,lay)
            fac010 = (1.  - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1.  - fs1) * fac01(iplon,lay)
            fac011 = (1.  - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng12
            tauself = selffac(iplon,lay)* (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac200 * absad(ind0+2,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig) + &
                    fac210 * absad(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                    (fac200 * absad(ind0-1,ig) + &
                    fac100 * absad(ind0,ig) + &
                    fac000 * absad(ind0+1,ig) + &
                    fac210 * absad(ind0+8,ig) + &
                    fac110 * absad(ind0+9,ig) + &
                    fac010 * absad(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac201 * absad(ind1+2,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig) + &
                    fac211 * absad(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                    (fac201 * absad(ind1-1,ig) + &
                    fac101 * absad(ind1,ig) + &
                    fac001 * absad(ind1+1,ig) + &
                    fac211 * absad(ind1+8,ig) + &
                    fac111 * absad(ind1+9,ig) + &
                    fac011 * absad(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig))
            endif

            taug(iplon,lay,ngs11+ig) = tau_major + tau_major1 &
                 + tauself + taufor
            fracsd(iplon,lay,ngs11+ig) = fracrefad(ig,jpl) + fpl * &
                 (fracrefad(ig,jpl+1)-fracrefad(ig,jpl))
         enddo
   
else
         do ig = 1, ng12
            taug(iplon,lay,ngs11+ig) = 0.0 
            fracsd(iplon,lay,ngs11+ig) = 0.0 
         enddo
    endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb12g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb13g( ncol, nlayers, taug, fracsd )
!----------------------------------------------------------------------------
!
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng13, ngs12
      use parrrtm, only : ngs12
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg13
! ------- Declarations -------
 real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
! Local 
      integer  :: lay, ind0, ind1, inds, indf, indm, ig
      integer  :: js, js1, jmco2, jmco, jpl
      real  :: speccomb, specparm, specmult, fs
      real  :: speccomb1, specparm1, specmult1, fs1
      real  :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real  :: speccomb_mco, specparm_mco, specmult_mco, fmco
      real  :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real  :: p, p4, fk0, fk1, fk2
      real  :: fac000, fac100, fac200, fac010, fac110, fac210
      real  :: fac001, fac101, fac201, fac011, fac111, fac211
      real  :: tauself, taufor, co2m1, co2m2, absco2 
      real  :: com1, com2, absco, abso3
      real  :: chi_co2, ratco2, adjfac, adjcolco2
      real  :: refrat_planck_a, refrat_m_a, refrat_m_a3
      real  :: tau_major, tau_major1
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
     iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Minor gas mapping levels :
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - co, p = 706 mb, t = 278.94 k
!     upper - o3, p = 95.5835 mb, t = 215.7 k

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

! P = 473.420 mb (Level 5)
      refrat_planck_a = chi_mlsd(1,5)/chi_mlsd(4,5)

! P = 1053. (Level 1)
      refrat_m_a = chi_mlsd(1,1)/chi_mlsd(4,1)

! P = 706. (Level 3)
      refrat_m_a3 = chi_mlsd(1,3)/chi_mlsd(4,3)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately.  

! Lower atmosphere loop
      if (lay <= laytrop(iplon)) then

         speccomb = colh2o(iplon,lay) + rat_h2on2o(iplon,lay)*coln2o(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(iplon,lay) + rat_h2on2o_1(iplon,lay)*coln2o(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mco2 = colh2o(iplon,lay) + refrat_m_a*coln2o(iplon,lay)
         specparm_mco2 = colh2o(iplon,lay)/speccomb_mco2
         if (specparm_mco2 .ge. oneminusd) specparm_mco2 = oneminusd
         specmult_mco2 = 8. *specparm_mco2
         jmco2 = 1 + int(specmult_mco2)
         fmco2 = mod(specmult_mco2,1.0 )

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(iplon,lay)/(coldry(iplon,lay))
         ratco2 = 1.e20 *chi_co2/3.55e-4 
         if (ratco2 .gt. 3.0 ) then
            adjfac = 2.0 +(ratco2-2.0 )**0.68 
            adjcolco2 = adjfac*3.55e-4*coldry(iplon,lay)*1.e-20 
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         speccomb_mco = colh2o(iplon,lay) + refrat_m_a3*coln2o(iplon,lay)
         specparm_mco = colh2o(iplon,lay)/speccomb_mco
         if (specparm_mco .ge. oneminusd) specparm_mco = oneminusd
         specmult_mco = 8. *specparm_mco
         jmco = 1 + int(specmult_mco)
         fmco = mod(specmult_mco,1.0 )

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*coln2o(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(13) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(13) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1.  - fs) * fac00(iplon,lay)
            fac010 = (1.  - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1.  - fs1) * fac01(iplon,lay)
            fac011 = (1.  - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng13
            tauself = selffac(iplon,lay)* (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor = forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            co2m1 = ka_mco2d(jmco2,indm,ig) + fmco2 * &
                 (ka_mco2d(jmco2+1,indm,ig) - ka_mco2d(jmco2,indm,ig))
            co2m2 = ka_mco2d(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco2d(jmco2+1,indm+1,ig) - ka_mco2d(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(iplon,lay) * (co2m2 - co2m1)
            com1 = ka_mcod(jmco,indm,ig) + fmco * &
                 (ka_mcod(jmco+1,indm,ig) - ka_mcod(jmco,indm,ig))
            com2 = ka_mcod(jmco,indm+1,ig) + fmco * &
                 (ka_mcod(jmco+1,indm+1,ig) - ka_mcod(jmco,indm+1,ig))
            absco = com1 + minorfrac(iplon,lay) * (com2 - com1)

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac200 * absad(ind0+2,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig) + &
                    fac210 * absad(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                    (fac200 * absad(ind0-1,ig) + &
                    fac100 * absad(ind0,ig) + &
                    fac000 * absad(ind0+1,ig) + &
                    fac210 * absad(ind0+8,ig) + &
                    fac110 * absad(ind0+9,ig) + &
                    fac010 * absad(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac201 * absad(ind1+2,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig) + &
                    fac211 * absad(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                    (fac201 * absad(ind1-1,ig) + &
                    fac101 * absad(ind1,ig) + &
                    fac001 * absad(ind1+1,ig) + &
                    fac211 * absad(ind1+8,ig) + &
                    fac111 * absad(ind1+9,ig) + &
                    fac011 * absad(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig))
            endif

            taug(iplon,lay,ngs12+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colco(iplon,lay)*absco
            fracsd(iplon,lay,ngs12+ig) = fracrefad(ig,jpl) + fpl * &
                 (fracrefad(ig,jpl+1)-fracrefad(ig,jpl))
         enddo
    else
         indm = indminor(iplon,lay)
         do ig = 1, ng13
            abso3 = kb_mo3d(indm,ig) + minorfrac(iplon,lay) * &
                 (kb_mo3d(indm+1,ig) - kb_mo3d(indm,ig))
            taug(iplon,lay,ngs12+ig) = colo3(iplon,lay)*abso3
            fracsd(iplon,lay,ngs12+ig) =  fracrefbd(ig)
         enddo
      endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb13g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb14g( ncol, nlayers , taug, fracsd)
!----------------------------------------------------------------------------
!
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng14, ngs13
      use parrrtm, only : ngs13
      use rrlw_kg14

! ------- Declarations -------
 real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
! Local 
      integer  :: lay, ind0, ind1, inds, indf, ig
      real  :: tauself, taufor
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
    iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum 
! and foreign continuum is interpolated (in temperature) separately.  

! Lower atmosphere loop
      if (lay <= laytrop(iplon)) then
         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(14) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(14) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         do ig = 1, ng14
            tauself = selffac(iplon,lay) * (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor =  forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            taug(iplon,lay,ngs13+ig) = colco2(iplon,lay) * &
                 (fac00(iplon,lay) * absad(ind0,ig) + &
                 fac10(iplon,lay) * absad(ind0+1,ig) + &
                 fac01(iplon,lay) * absad(ind1,ig) + &
                 fac11(iplon,lay) * absad(ind1+1,ig)) &
                 + tauself + taufor
            fracsd(iplon,lay,ngs13+ig) = fracrefad(ig)
         enddo
    else
         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(14) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(14) + 1
         do ig = 1, ng14
            taug(iplon,lay,ngs13+ig) = colco2(iplon,lay) * &
                 (fac00(iplon,lay) * absbd(ind0,ig) + &
                 fac10(iplon,lay) * absbd(ind0+1,ig) + &
                 fac01(iplon,lay) * absbd(ind1,ig) + &
                 fac11(iplon,lay) * absbd(ind1+1,ig))
            fracsd(iplon,lay,ngs13+ig) = fracrefbd(ig)
         enddo
      endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb14g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb15g( ncol, nlayers , taug, fracsd)
!----------------------------------------------------------------------------
!
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng15, ngs14
      use parrrtm, only : ngs14
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg15

! ------- Declarations -------
 real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
! Local 
      integer  :: lay, ind0, ind1, inds, indf, indm, ig
      integer  :: js, js1, jmn2, jpl
      real  :: speccomb, specparm, specmult, fs
      real  :: speccomb1, specparm1, specmult1, fs1
      real  :: speccomb_mn2, specparm_mn2, specmult_mn2, fmn2
      real  :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real  :: p, p4, fk0, fk1, fk2
      real  :: fac000, fac100, fac200, fac010, fac110, fac210
      real  :: fac001, fac101, fac201, fac011, fac111, fac211
      real  :: scalen2, tauself, taufor, n2m1, n2m2, taun2 
      real  :: refrat_planck_a, refrat_m_a
      real  :: tau_major, tau_major1
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
       iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif
! Minor gas mapping level : 
!     Lower - Nitrogen Continuum, P = 1053., T = 294.

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower atmosphere.
! P = 1053. mb (Level 1)
      refrat_planck_a = chi_mlsd(4,1)/chi_mlsd(2,1)

! P = 1053.
      refrat_m_a = chi_mlsd(4,1)/chi_mlsd(2,1)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately.  

! Lower atmosphere loop
   if (lay <= laytrop(iplon)) then

         speccomb = coln2o(iplon,lay) + rat_n2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = coln2o(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = coln2o(iplon,lay) + rat_n2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = coln2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_mn2 = coln2o(iplon,lay) + refrat_m_a*colco2(iplon,lay)
         specparm_mn2 = coln2o(iplon,lay)/speccomb_mn2
         if (specparm_mn2 .ge. oneminusd) specparm_mn2 = oneminusd
         specmult_mn2 = 8. *specparm_mn2
         jmn2 = 1 + int(specmult_mn2)
         fmn2 = mod(specmult_mn2,1.0 )

         speccomb_planck = coln2o(iplon,lay)+refrat_planck_a*colco2(iplon,lay)
         specparm_planck = coln2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(15) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(15) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)
         
         scalen2 = colbrd(iplon,lay)*scaleminor(iplon,lay)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1.  - fs) * fac00(iplon,lay)
            fac010 = (1.  - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif
         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1.  - fs1) * fac01(iplon,lay)
            fac011 = (1.  - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng15
            tauself = selffac(iplon,lay)* (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor =  forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 
            n2m1 = ka_mn2d(jmn2,indm,ig) + fmn2 * &
                 (ka_mn2d(jmn2+1,indm,ig) - ka_mn2d(jmn2,indm,ig))
            n2m2 = ka_mn2d(jmn2,indm+1,ig) + fmn2 * &
                 (ka_mn2d(jmn2+1,indm+1,ig) - ka_mn2d(jmn2,indm+1,ig))
            taun2 = scalen2 * (n2m1 + minorfrac(iplon,lay) * (n2m2 - n2m1))

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac200 * absad(ind0+2,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig) + &
                    fac210 * absad(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                    (fac200 * absad(ind0-1,ig) + &
                    fac100 * absad(ind0,ig) + &
                    fac000 * absad(ind0+1,ig) + &
                    fac210 * absad(ind0+8,ig) + &
                    fac110 * absad(ind0+9,ig) + &
                    fac010 * absad(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig))
            endif 

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac201 * absad(ind1+2,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig) + &
                    fac211 * absad(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                    (fac201 * absad(ind1-1,ig) + &
                    fac101 * absad(ind1,ig) + &
                    fac001 * absad(ind1+1,ig) + &
                    fac211 * absad(ind1+8,ig) + &
                    fac111 * absad(ind1+9,ig) + &
                    fac011 * absad(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig))
            endif

            taug(iplon,lay,ngs14+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + taun2
            fracsd(iplon,lay,ngs14+ig) = fracrefad(ig,jpl) + fpl * &
                 (fracrefad(ig,jpl+1)-fracrefad(ig,jpl))
         enddo
    
    else
         do ig = 1, ng15
            taug(iplon,lay,ngs14+ig) = 0.0 
            fracsd(iplon,lay,ngs14+ig) = 0.0 
         enddo
      endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb15g

!----------------------------------------------------------------------------
      _gpuker subroutine taugb16g( ncol, nlayers , taug, fracsd)
!----------------------------------------------------------------------------
!
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!----------------------------------------------------------------------------

! ------- Modules -------

!      use parrrtm, only : ng16, ngs15
      use parrrtm, only : ngs15
      use rrlw_ref, only : chi_mlsd
      use rrlw_kg16

! ------- Declarations -------
 real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
     
! Local 
      integer  :: lay, ind0, ind1, inds, indf, ig
      integer  :: js, js1, jpl
      real  :: speccomb, specparm, specmult, fs
      real  :: speccomb1, specparm1, specmult1, fs1
      real  :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real  :: p, p4, fk0, fk1, fk2
      real  :: fac000, fac100, fac200, fac010, fac110, fac210
      real  :: fac001, fac101, fac201, fac011, fac111, fac211
      real  :: tauself, taufor
      real  :: refrat_planck_a
      real  :: tau_major, tau_major1
      integer , value, intent(in) :: ncol, nlayers
      integer  :: iplon
#ifdef _CUDA
     iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      lay = (blockidx%y-1) * blockdim%y + threadidx%y
      if (iplon <= ncol .and. lay <= nlayers) then
#else
      do iplon = 1, ncol
      do lay = 1, nlayers
#endif 
! Calculate reference ratio to be used in calculation of Planck
! fraction in lower atmosphere.

! P = 387. mb (Level 6)
      refrat_planck_a = chi_mlsd(1,6)/chi_mlsd(6,6)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature,and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately.  

! Lower atmosphere loop
 if (lay <= laytrop(iplon)) then
         speccomb = colh2o(iplon,lay) + rat_h2och4(iplon,lay)*colch4(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminusd) specparm = oneminusd
         specmult = 8. *(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0 )

         speccomb1 = colh2o(iplon,lay) + rat_h2och4_1(iplon,lay)*colch4(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminusd) specparm1 = oneminusd
         specmult1 = 8. *(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0 )

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colch4(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminusd) specparm_planck=oneminusd
         specmult_planck = 8. *specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0 )

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspad(16) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspad(16) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)

         if (specparm .lt. 0.125 ) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875 ) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1.  - fs) * fac00(iplon,lay)
            fac010 = (1.  - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125 ) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875 ) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0 *p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1.  - fs1) * fac01(iplon,lay)
            fac011 = (1.  - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng16
            tauself = selffac(iplon,lay)* (selfrefd(inds,ig) + selffrac(iplon,lay) * &
                 (selfrefd(inds+1,ig) - selfrefd(inds,ig)))
            taufor =  forfac(iplon,lay) * (forrefd(indf,ig) + forfrac(iplon,lay) * &
                 (forrefd(indf+1,ig) - forrefd(indf,ig))) 

            if (specparm .lt. 0.125 ) then
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac200 * absad(ind0+2,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig) + &
                    fac210 * absad(ind0+11,ig))
            else if (specparm .gt. 0.875 ) then
               tau_major = speccomb * &
                    (fac200 * absad(ind0-1,ig) + &
                    fac100 * absad(ind0,ig) + &
                    fac000 * absad(ind0+1,ig) + &
                    fac210 * absad(ind0+8,ig) + &
                    fac110 * absad(ind0+9,ig) + &
                    fac010 * absad(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absad(ind0,ig) + &
                    fac100 * absad(ind0+1,ig) + &
                    fac010 * absad(ind0+9,ig) + &
                    fac110 * absad(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125 ) then
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac201 * absad(ind1+2,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig) + &
                    fac211 * absad(ind1+11,ig))
            else if (specparm1 .gt. 0.875 ) then
               tau_major1 = speccomb1 * &
                    (fac201 * absad(ind1-1,ig) + &
                    fac101 * absad(ind1,ig) + &
                    fac001 * absad(ind1+1,ig) + &
                    fac211 * absad(ind1+8,ig) + &
                    fac111 * absad(ind1+9,ig) + &
                    fac011 * absad(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absad(ind1,ig) + &
                    fac101 * absad(ind1+1,ig) + &
                    fac011 * absad(ind1+9,ig) + &
                    fac111 * absad(ind1+10,ig))
            endif

            taug(iplon,lay,ngs15+ig) = tau_major + tau_major1 &
                 + tauself + taufor
            fracsd(iplon,lay,ngs15+ig) = fracrefad(ig,jpl) + fpl * &
                 (fracrefad(ig,jpl+1)-fracrefad(ig,jpl))
         enddo
else
         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspbd(16) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspbd(16) + 1
         do ig = 1, ng16
            taug(iplon,lay,ngs15+ig) = colch4(iplon,lay) * &
                 (fac00(iplon,lay) * absbd(ind0,ig) + &
                 fac10(iplon,lay) * absbd(ind0+1,ig) + &
                 fac01(iplon,lay) * absbd(ind1,ig) + &
                 fac11(iplon,lay) * absbd(ind1+1,ig))
            fracsd(iplon,lay,ngs15+ig) = fracrefbd(ig)
         enddo
     endif

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif

      end subroutine taugb16g

      _gpuker subroutine addAerosols( ncol, nlayers, ngptlw, nbndlw, ngbd, taug)

      integer , intent(in), value :: ncol, nlayers, ngptlw, nbndlw
        integer , intent(in) :: ngbd(:)
        
        
        integer  :: iplon, lay, ig
         real  _gpudev :: taug(:,:,:)
     
#ifdef _CUDA     
        iplon = (blockidx%x-1) * blockdim%x + threadidx%x
        lay = (blockidx%y-1) * blockdim%y + threadidx%y
        ig = (blockidx%z-1) * blockdim%z + threadidx%z
        if (iplon<=ncol .and. lay<=nlayers .and. ig<=ngptlw) then
#else
        do iplon = 1, ncol
        do lay = 1, nlayers
        do ig = 1, ngptlw
#endif



            taug(iplon, lay, ig) = taug(iplon, lay, ig) + tauaa(iplon, lay, ngbd(ig))

#ifdef _CUDA
        endif
#else
     end do
     end do
     end do
#endif
        




      end subroutine

!----------------------------------------------------------------------------
      subroutine taumolg(iplon, ncol, nlayers, ngbd, taug, fracsd)
!----------------------------------------------------------------------------

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
! *     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
! *     COMMON /PLANKG/   fracsd(MXLAY,MG)                                       *
! *                                                                             *
! *  Input                                                                      *
! *                                                                             *
! *     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
! *     COMMON /PRECISE/  oneminusd                                              *
! *     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
! *     &                 PZ(0:MXLAY),TZ(0:MXLAY)                               *
! *     COMMON /PROFDATA/ LAYTROP,                                              *
! *    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),             *
! *    &                  COLN2O(MXLAY),colco(MXLAY),COLCH4(MXLAY),             *
! *    &                  COLO2(MXLAY)
! *     COMMON /INTFAC/   fac00(iplon,MXLAY),fac01(iplon,MXLAY),                            *
! *    &                  FAC10(MXLAY),fac11(iplon,MXLAY)                             *
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
! *     oneminusd - since problems are caused in some cases by interpolation     *
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
!*******************************************************************************

      use parrrtm, only : ng1


! ------- Declarations -------

! ----- Input -----
      integer , intent(in) :: iplon           ! the column number (move to calculated in kernel)
      integer , intent(in) :: ncol            ! the total number of columns
      integer , intent(in) :: nlayers         ! total number of layers
      integer  _gpudev, intent(in) :: ngbd(:)
      real , intent(in) _gpudev :: fracsd(:,:,:)
      real , intent(in) _gpudev :: taug(:,:,:)
   
      !real  :: taugcc(ncol, nlayers, 140)

! ----- Output -----
  
      integer :: i,j,err
      real :: t1, t2
#ifdef _CUDA
      type(dim3) :: dimGrid, dimBlock
    
#endif
#ifdef _CUDA
      !dimGrid = dim3( (ncol + 127) / 128, 1, 1)
	  !dimBlock = dim3( 128,1,1)

      dimGrid = dim3( (ncol + 63) / 64, ((nlayers+1)/2), 1)
      dimBlock = dim3( 64, 2, 1)
      
#endif   


! Calculate gaseous optical depth and planck fractions for each spectral band.

! (dmb 2012) Here we configure the grid and thread blocks.  These subroutines are
! only parallelized across the column dimension so the blocks are one dimensional.

! (dmb 2012) Call all 16 kernels in sequence



      call taugb1g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb2g _gpuchv (ncol, nlayers, taug, fracsd)

      
      call taugb3g _gpuchv (ncol, nlayers, taug, fracsd)

  
      call taugb4g _gpuchv (ncol, nlayers, taug, fracsd)

      
      call taugb5g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb6g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb7g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb8g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb9g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb10g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb11g _gpuchv (ncol, nlayers, taug, fracsd)
 
      call taugb12g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb13g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb14g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb15g _gpuchv (ncol, nlayers, taug, fracsd)

      call taugb16g _gpuchv (ncol, nlayers, taug, fracsd)
      
     

#ifdef _CUDA
      dimGrid = dim3( (ncol+ 255) / 256, nlayers, ngptlw )
      dimBlock = dim3( 256, 1, 1)
#endif

      ! (dmb 2012) This code used to be in the main rrtmg_lw_rad source file
      ! We add the aerosol optical depths to the gas optical depths
      call addAerosols _gpuchv (ncol, nlayers, ngptlw, nbndlw, ngbd, taug)

      end subroutine taumolg


      ! (dmb 2012) Allocate all of the needed memory for the taumol subroutines
      subroutine allocateGPUTaumol(ncol, nlayers, npart)

          integer , intent(in) :: ncol
          integer , intent(in) :: nlayers
          integer , intent(in) :: npart
          integer :: i

           sreg( wx1 , ncol, nlayers )
          sreg( wx2 , ncol, nlayers )
          sreg( wx3 , ncol, nlayers )
          sreg( wx4 , ncol, nlayers )

         
         
          sreg( jp , ncol, nlayers )
          sreg( jt , ncol, nlayers )
          sreg( jt1 , ncol, nlayers )
          sreg( colh2o , ncol, nlayers )
          sreg( colco2 , ncol, nlayers )
          sreg( colo3 , ncol, nlayers )
          sreg( coln2o , ncol, nlayers )
          sreg( colco , ncol, nlayers )
          sreg( colch4 , ncol, nlayers )
          sreg( colo2 , ncol, nlayers )
          sreg( colbrd , ncol, nlayers )
          sreg( indself , ncol, nlayers )
          sreg( indfor , ncol, nlayers )
          sreg( selffac , ncol, nlayers )
          sreg( selffrac , ncol, nlayers )
          sreg( forfac , ncol, nlayers )
          sreg( forfrac , ncol, nlayers )
          sreg( indminor , ncol, nlayers )
          sreg( minorfrac , ncol, nlayers )
          sreg( scaleminor , ncol, nlayers )
          sreg( scaleminorn2 , ncol, nlayers )
        
          sreg( fac00 , ncol, nlayers )
          sreg( fac10 , ncol, nlayers )
          sreg( fac01 , ncol, nlayers )
          sreg( fac11 , ncol, nlayers )
          sreg( rat_h2oco2 , ncol, nlayers )
          sreg( rat_h2oco2_1 , ncol, nlayers )
          sreg( rat_h2oo3 , ncol, nlayers )
          sreg( rat_h2oo3_1 , ncol, nlayers )
          sreg( rat_h2on2o , ncol, nlayers )
          sreg( rat_h2on2o_1 , ncol, nlayers )
          sreg( rat_h2och4 , ncol, nlayers )
          sreg( rat_h2och4_1 , ncol, nlayers )
          sreg( rat_n2oco2 , ncol, nlayers )
          sreg( rat_n2oco2_1 , ncol, nlayers )
          sreg( rat_o3co2 , ncol, nlayers )
          sreg( rat_o3co2_1 , ncol, nlayers )
#ifdef _CUDA
            call dflush()
#endif
          allocate( pavel( ncol, nlayers ))
          dreg( wx1 , ncol, nlayers )
          dreg( wx2 , ncol, nlayers )
          dreg( wx3 , ncol, nlayers )
          dreg( wx4 , ncol, nlayers )

          allocate( coldry( ncol, nlayers ))
         
          dreg( jp , ncol, nlayers )
          dreg( jt , ncol, nlayers )
          dreg( jt1 , ncol, nlayers )
          dreg( colh2o , ncol, nlayers )
          dreg( colco2 , ncol, nlayers )
          dreg( colo3 , ncol, nlayers )
          dreg( coln2o , ncol, nlayers )
          dreg( colco , ncol, nlayers )
          dreg( colch4 , ncol, nlayers )
          dreg( colo2 , ncol, nlayers )
          dreg( colbrd , ncol, nlayers )
          dreg( indself , ncol, nlayers )
          dreg( indfor , ncol, nlayers )
          dreg( selffac , ncol, nlayers )
          dreg( selffrac , ncol, nlayers )
          dreg( forfac , ncol, nlayers )
          dreg( forfrac , ncol, nlayers )
          dreg( indminor , ncol, nlayers )
          dreg( minorfrac , ncol, nlayers )
          dreg( scaleminor , ncol, nlayers )
          dreg( scaleminorn2 , ncol, nlayers )

          dreg( fac00 , ncol, nlayers )
          dreg( fac10 , ncol, nlayers )
          dreg( fac01 , ncol, nlayers )
          dreg( fac11 , ncol, nlayers )
          dreg( rat_h2oco2 , ncol, nlayers )
          dreg( rat_h2oco2_1 , ncol, nlayers )
          dreg( rat_h2oo3 , ncol, nlayers )
          dreg( rat_h2oo3_1 , ncol, nlayers )
          dreg( rat_h2on2o , ncol, nlayers )
          dreg( rat_h2on2o_1 , ncol, nlayers )
          dreg( rat_h2och4 , ncol, nlayers )
          dreg( rat_h2och4_1 , ncol, nlayers )
          dreg( rat_n2oco2 , ncol, nlayers )
          dreg( rat_n2oco2_1 , ncol, nlayers )
          dreg( rat_o3co2 , ncol, nlayers )
          dreg( rat_o3co2_1 , ncol, nlayers )


          allocate( laytrop( ncol ))
          allocate( tauaa( ncol, nlayers, nbndlw ))
          allocate( nspad( nbndlw ))
          allocate( nspbd( nbndlw ))
        
      end subroutine

      ! (dmb 2012) Perform the necessary cleanup of the GPU arrays
      subroutine deallocateGPUTaumol()

#ifdef _CUDA
          call dbclean
          call dclean
#endif
          deallocate( pavel)
      
          deallocate( tauaa )
          deallocate( laytrop)
       
          deallocate( nspad)
          deallocate( nspbd)
          deallocate( coldry)
      end subroutine

       
      subroutine copyGPUTaumolMol( colstart, pncol, nlayers, colh2oc, colco2c, colo3c, coln2oc, colch4c, colo2c,&
                                   px1,px2,px3,px4, npart)
        
        integer, value, intent(in) :: colstart, pncol, nlayers, npart
        real , intent(in) :: colh2oc(:,:), colco2c(:,:), colo3c(:,:), coln2oc(:,:), &
                                     colch4c(:,:), colo2c(:,:), px1(:,:), px2(:,:), px3(:,:), px4(:,:)

        if (npart > 1) then
        colh2o(1:pncol, :) = colh2oc( colstart:(colstart+pncol-1), 1:nlayers)
        colco2(1:pncol, :) = colco2c( colstart:(colstart+pncol-1), 1:nlayers)
        colo3(1:pncol, :) = colo3c( colstart:(colstart+pncol-1), 1:nlayers)
        coln2o(1:pncol, :) = coln2oc( colstart:(colstart+pncol-1), 1:nlayers)
      
        colch4(1:pncol, :) = colch4c( colstart:(colstart+pncol-1), 1:nlayers)
        colo2(1:pncol, :) = colo2c( colstart:(colstart+pncol-1), 1:nlayers)
        wx1(1:pncol, :) = px1(colstart:(colstart+pncol-1), 1:nlayers)
        wx2(1:pncol, :) = px2(colstart:(colstart+pncol-1), 1:nlayers)
        wx3(1:pncol, :) = px3(colstart:(colstart+pncol-1), 1:nlayers)
        wx4(1:pncol, :) = px4(colstart:(colstart+pncol-1), 1:nlayers)
        else
         colh2o = colh2oc
         colco2 = colco2c
         colo3 = colo3c
         coln2o = coln2oc
         colch4 = colch4c
         colo2 = colo2c
         wx1 = px1
         wx2 = px2
         wx3 = px3
         wx4 = px4

        endif
        colco = 0
      end subroutine

      ! (dmb 2012) Copy the needed data from the CPU to the GPU.  I had to separate this
      ! out into 16 separate functions to correspond with the 16 taumol subroutines.
      subroutine copyGPUTaumol(pavelc, coldryc, tauap, pncol, colstart, nlay, npart)
      use rrlw_kg01, only : copyToGPU1, reg1
      use rrlw_kg02, only : copyToGPU2, reg2
      use rrlw_kg03, only : copyToGPU3, reg3
      use rrlw_kg04, only : copyToGPU4, reg4
      use rrlw_kg05, only : copyToGPU5, reg5
      use rrlw_kg06, only : copyToGPU6, reg6
      use rrlw_kg07, only : copyToGPU7, reg7
      use rrlw_kg08, only : copyToGPU8, reg8
      use rrlw_kg09, only : copyToGPU9, reg9
      use rrlw_kg10, only : copyToGPU10, reg10
      use rrlw_kg11, only : copyToGPU11, reg11
      use rrlw_kg12, only : copyToGPU12, reg12
      use rrlw_kg13, only : copyToGPU13, reg13
      use rrlw_kg14, only : copyToGPU14, reg14
      use rrlw_kg15, only : copyToGPU15, reg15
      use rrlw_kg16, only : copyToGPU16, reg16
      use rrlw_ref, only  : copyToGPUref
      real , intent(in) :: pavelc(:,:)           ! layer pressures (mb) 
                                                      !    Dimensions: (nlayers)
      !real  , intent(in) :: wxc(:,:,:)            ! cross-section amounts (mol/cm2)
                                                      !    Dimensions: (maxxsec,nlayers)
      real  , intent(in) :: coldryc(:,:)          ! column amount (dry air)
                                                      !    Dimensions: (nlayers)

      real , intent(in) :: tauap(:,:,:)
      integer, intent(in)      :: pncol, colstart, nlay, npart
     
      
    
      call reg1
      call reg2
      call reg3
      call reg4
      call reg5
      call reg6
      call reg7
      call reg8
      call reg9
      call reg10
      call reg11
      call reg12
      call reg13
      call reg14
      call reg15
      call reg16
      
      dbflushreg()
      call CopyToGPU1
      call CopyToGPU2
      call CopyToGPU3
      call CopyToGPU4
      call CopyToGPU5
      call CopyToGPU6
      call CopyToGPU7
      call CopyToGPU8
      call CopyToGPU9
      call CopyToGPU10
      call CopyToGPU11
      call CopyToGPU12
      call CopyToGPU13
      call CopyToGPU14
      call CopyToGPU15
      call CopyToGPU16

      nspad= nspa
      nspbd= nspb
      pavel= pavelc
      coldry= coldryc
      

      oneminusd = oneminus

      dbflushcop()
     
      if (npart > 1) then
      tauaa(1:pncol, :, :)  = tauap(colstart:(colstart+pncol-1), :, :)
      else
      tauaa = tauap
      endif
 
      end subroutine 

      end module gpu_rrtmg_lw_taumol

