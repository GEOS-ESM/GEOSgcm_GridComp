
!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

#ifdef _CUDA
#define gpu_device ,device
#else
#define gpu_device 
#endif

      module rrtmg_sw_spcvmc

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
      use parrrsw, only : nbndsw, ngptsw, mxmol, jpband, mxlay
      use rrsw_tbl, only : tblint, bpade, od_lo, exp_tbl
      use rrsw_vsn, only : hvrspc, hnamspc
      use rrsw_wvn, only : ngc, ngs, ngb
      
      use rrtmg_sw_taumol, only: taumol_sw
     

      implicit none

      contains

! ---------------------------------------------------------------------------
      subroutine spcvmc_sw &
            (cc,tncol, ncol, nlayers, istart, iend, icpr, idelm, iout, &
             pavel, tavel, pz, tz, tbound, palbd, palbp, &
             pcldfmc, ptaucmc, pasycmc, pomgcmc, ptaormc, &
             ptaua, pasya, pomga, prmu0, coldry,  adjflux, &
             isolvar, svar_f, svar_s, svar_i, &
             svar_f_bnd, svar_s_bnd, svar_i_bnd, &
             laytrop, layswtch, laylow, jp, jt, jt1, &
             co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
             fac00, fac01, fac10, fac11, &
             selffac, selffrac, indself, forfac, forfrac, indfor, &
             pbbfd, pbbfu, pbbcd, pbbcu, puvfd, puvcd, pnifd, pnicd, &
             pbbfddir, pbbcddir, puvfddir, puvcddir, pnifddir, pnicddir,&
            zgco,zomco,zrdnd,zref,zrefo,zrefd,zrefdo,ztauo,zdbt,ztdbt,&
             ztra,ztrao,ztrad,ztrado,zfd,zfu,ztaug, ztaur, zsflxzen, ssi,&
           znirr,znirf,zparr,zparf,zuvrr,zuvrf)
! ---------------------------------------------------------------------------
!
! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
!          using the two-stream method of H. Barker and McICA, the Monte-Carlo
!          Independent Column Approximation, for the representation of 
!          sub-grid cloud variability (i.e. cloud overlap).
!
! Interface:  *spcvmc_sw* is called from *rrtmg_sw.F90* or rrtmg_sw.1col.F90*
!
! Method:
!    Adapted from two-stream model of H. Barker;
!    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
!        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
!
! Modifications:
!
! Original: H. Barker
! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
! Revision: Code modified so that delta scaling is not done in cloudy profiles
!           if routine cldprop is used; delta scaling can be applied by swithcing
!           code below if cldprop is not used to get cloud properties. 
!           AER, Jan 2005
! Revision: Modified to use McICA: MJIacono, AER, Nov 2005
! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
!           Aug 2007 
!
! ------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------



      integer , intent(in) :: tncol, ncol,cc
      integer, intent(in) :: nlayers
      integer, intent(in) :: istart
      integer, intent(in) :: iend
      integer, intent(in) :: icpr
      integer, intent(in) :: idelm   ! delta-m scaling flag
                                              ! [0 = direct and diffuse fluxes are unscaled]
                                              ! [1 = direct and diffuse fluxes are scaled]
      integer, intent(in) :: iout
      integer , intent(in) :: laytrop(:)
      integer , intent(in) :: layswtch(:)
      integer , intent(in) :: laylow(:)

      integer , intent(in) :: indfor(:,:) 
                                                               !   Dimensions: (nlayers)
      integer , intent(in) :: indself(:,:) 
                                                               !   Dimensions: (nlayers)
      integer , intent(in) :: jp(:,:) 
                                                               !   Dimensions: (nlayers)
      integer , intent(in) :: jt(:,:) 
                                                               !   Dimensions: (nlayers)
      integer , intent(in) :: jt1(:,:) 
                                                               !   Dimensions: (nlayers)

      real , intent(in) :: pavel(:,:)                     ! layer pressure (hPa, mb) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: tavel(:,:)                     ! layer temperature (K)
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: pz(:,0:)                       ! level (interface) pressure (hPa, mb)
                                                               !   Dimensions: (0:nlayers)
      real , intent(in) :: tz(:,0:)                       ! level temperatures (hPa, mb)
                                                               !   Dimensions: (0:nlayers)
      real , intent(in) :: tbound(:)                      ! surface temperature (K)
      real , intent(in) :: coldry(:,:)                    ! dry air column density (mol/cm2)
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: colmol(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: adjflux(:)                  ! Earth/Sun distance adjustment
                                                               !   Dimensions: (jpband)
! Solar variability
      integer, intent(in) :: isolvar                  ! Flag for solar variability method
      real, intent(in) :: svar_f                      ! Solar variability facular multiplier
      real, intent(in) :: svar_s                      ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                      ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(jpband)          ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(jpband)          ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(jpband)          ! Solar variability baseline irradiance multiplier (by band)

      real , intent(in) :: palbd(:,:)                     ! surface albedo (diffuse)
                                                               !   Dimensions: (nbndsw)
      real , intent(in) :: palbp(:,:)                     ! surface albedo (direct)
                                                               !   Dimensions: (nbndsw)
      real , intent(in) :: prmu0(:)                       ! cosine of solar zenith angle
      real , intent(in) :: pcldfmc(:,:,:)                 ! cloud fraction [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real , intent(in) :: ptaucmc(:,:,:)                 ! cloud optical depth [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real , intent(in) :: pasycmc(:,:,:)                 ! cloud asymmetry parameter [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real , intent(in) :: pomgcmc(:,:,:)                 ! cloud single scattering albedo [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real , intent(in) :: ptaormc(:,:,:)                 ! cloud optical depth, non-delta scaled [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
   
                                                               !   Dimensions: (nlayers,ngptsw)
      real , intent(in) :: ptaua(:,:,:)                  ! aerosol optical depth
                                                               !   Dimensions: (nlayers,nbndsw)
      real , intent(in) :: pasya(:,:,:)                  ! aerosol asymmetry parameter
                                                               !   Dimensions: (nlayers,nbndsw)
      real , intent(in) :: pomga(:,:,:)                  ! aerosol single scattering albedo
                                                               !   Dimensions: (nlayers,nbndsw)

                                                               
                                                               
      real , intent(in) :: colh2o(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: colco2(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: colch4(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: co2mult(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: colo3(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: colo2(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: coln2o(:,:) 
                                                               !   Dimensions: (nlayers)

      real , intent(in) :: forfac(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: forfrac(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: selffac(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: selffrac(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: fac00(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: fac01(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: fac10(:,:) 
                                                               !   Dimensions: (nlayers)
      real , intent(in) :: fac11(:,:) 
                                                               !   Dimensions: (nlayers)

      real, intent(inout) gpu_device :: zgco(tncol,ngptsw,nlayers+1), zomco(tncol,ngptsw,nlayers+1)  
      real, intent(inout) gpu_device  :: zrdnd(tncol,ngptsw,nlayers+1) 
      real, intent(inout) gpu_device  :: zref(tncol,ngptsw,nlayers+1)  , zrefo(tncol,ngptsw,nlayers+1)  
      real, intent(inout) gpu_device  :: zrefd(tncol,ngptsw,nlayers+1)  , zrefdo(tncol,ngptsw,nlayers+1)  
      real, intent(inout) gpu_device  :: ztauo(tncol,ngptsw,nlayers)  
      real, intent(inout) gpu_device  :: zdbt(tncol,ngptsw,nlayers+1)  ,ztdbt(tncol,ngptsw,nlayers+1)   
      real, intent(inout) gpu_device  :: ztra(tncol,ngptsw,nlayers+1)  , ztrao(tncol,ngptsw,nlayers+1)  
      real, intent(inout) gpu_device  :: ztrad(tncol,ngptsw,nlayers+1)  , ztrado(tncol,ngptsw,nlayers+1)  
      real, intent(inout) gpu_device  :: zfd(tncol,ngptsw,nlayers+1)  , zfu(tncol,ngptsw,nlayers+1)   
      real, intent(inout) gpu_device :: ztaur(tncol,nlayers,ngptsw), ztaug(tncol,nlayers,ngptsw) 
      real, intent(inout) gpu_device :: zsflxzen(tncol,ngptsw)
      real, intent(inout) gpu_device :: ssi(tncol,ngptsw)
   

! ------- Output -------
                                                               !   All Dimensions: (nlayers+1)
      real , intent(out) :: pbbcd(:,:) 
      real , intent(out) :: pbbcu(:,:) 
      real , intent(out) :: pbbfd(:,:) 
      real , intent(out) :: pbbfu(:,:) 
      real , intent(out) :: pbbfddir(:,:) 
      real , intent(out) :: pbbcddir(:,:) 

      real , intent(out) :: puvcd(:,:) 
      real , intent(out) :: puvfd(:,:) 
      real , intent(out) :: puvcddir(:,:) 
      real , intent(out) :: puvfddir(:,:) 

      real , intent(out) :: pnicd(:,:) 
      real , intent(out) :: pnifd(:,:) 
      real , intent(out) :: pnicddir(:,:) 
      real , intent(out) :: pnifddir(:,:) 
      
      real, intent(out), dimension(:) :: znirr,znirf,zparr,zparf,zuvrr,zuvrf


! ------- Local -------


      integer  :: klev
      integer  :: ibm, ikl, ikp, ikx
      integer :: iw, jb, jg, jl, jk

      integer :: itind

      real :: tblind, ze1
      real :: zclear, zcloud

        
      real  :: zincflx, ze2
     

      real :: zdbtmc, zdbtmo, zf, zgw, zreflect
      real :: zwf, tauorig, repclc



! Arrays from rrtmg_sw_vrtqdr routine


  
      integer :: iplon


! ------------------------------------------------------------------




 !$acc kernels     
         pbbcd =0. 
         pbbcu =0. 
         pbbfd =0. 
         pbbfu =0. 
         pbbcddir =0. 
         pbbfddir =0. 
         puvcd =0. 
         puvfd =0. 
         puvcddir =0. 
         puvfddir =0. 
         pnicd =0. 
         pnifd =0. 
         pnicddir =0. 
         pnifddir =0.
         zsflxzen = 0.
         ssi = 0.
         znirr=0.
         znirf=0.
         zparr=0.
         zparf=0.
         zuvrr=0.
         zuvrf=0.
      klev = nlayers
!$acc end kernels      






         
! Calculate the optical depths for gaseous absorption and Rayleigh scattering     
      call taumol_sw(ncol,nlayers, &
                     colh2o , colco2 , colch4 , colo2 , &
                     colo3 , colmol , &
                     laytrop, jp, jt, jt1, &
                     fac00, fac01, fac10, fac11, &
                     selffac , selffrac , indself , forfac , forfrac ,&
                     indfor , &
                     isolvar, svar_f, svar_s, svar_i, &
                     svar_f_bnd, svar_s_bnd, svar_i_bnd, &
                     ssi, zsflxzen , ztaug, ztaur)



   repclc = 1.e-12 


!$acc kernels 
do iplon = 1, ncol

   
! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

    do iw = 1, 112
          jb = ngb(iw)
          ibm = jb-15


! Clear-sky    
!   TOA direct beam    
           

            
           
! Cloudy-sky    
!   Surface values
            ztrao(iplon,iw,klev+1)   =0.0 
            ztrado(iplon,iw,klev+1)  =0.0 
            zrefo(iplon,iw,klev+1)   =palbp(iplon,ibm) 
            zrefdo(iplon,iw,klev+1)  =palbd(iplon,ibm) 
           
! Total sky    
!   TOA direct beam    
            ztdbt(iplon,iw,1)  =1.0 
    
!   Surface values
            zdbt(iplon,iw,klev+1)   =0.0 
            ztra(iplon,iw,klev+1)   =0.0 
            ztrad(iplon,iw,klev+1)  =0.0 
            zref(iplon,iw,klev+1)   =palbp(iplon,ibm) 
            zrefd(iplon,iw,klev+1)  =palbd(iplon,ibm) 
    

    end do
end do
!$acc end kernels     



!$acc kernels loop 
do iplon = 1, ncol
    !$acc loop private(zf, zwf, ibm, ikl, jb)
    do iw = 1, 112
            !$acc loop seq
            do jk=1,klev

               ikl=klev+1-jk
                jb = ngb(iw)
               ibm = jb-15

! Clear-sky optical parameters including aerosols
               ztauo(iplon,iw,jk)   = ztaur(iplon,ikl,iw)  + ztaug(iplon,ikl,iw) + ptaua(iplon,ikl,ibm)      
               zomco(iplon,iw,jk)   = ztaur(iplon,ikl,iw) + ptaua(iplon,ikl,ibm) * pomga(iplon,ikl,ibm)
               zgco(iplon,iw,jk) = pasya(iplon,ikl,ibm) * pomga(iplon,ikl,ibm) * ptaua(iplon,ikl,ibm) / zomco(iplon,iw,jk)   
               zomco(iplon,iw,jk)   = zomco(iplon,iw,jk) / ztauo(iplon,iw,jk)
               
               zf = zgco(iplon, iw, jk)
               zf = zf * zf
               zwf = zomco(iplon, iw, jk) * zf
               ztauo(iplon, iw, jk) = (1.0 - zwf) * ztauo(iplon, iw, jk)
               zomco(iplon, iw, jk) = (zomco(iplon, iw, jk) - zwf) / (1.0 - zwf)
               zgco(iplon, iw, jk) = (zgco(iplon, iw, jk) - zf) / (1.0 - zf)



            enddo    
      end do
end do
!$acc end kernels               


! Clear sky reflectivities
            call reftra_sw (ncol, nlayers, &
                            pcldfmc, zgco, prmu0, ztauo, zomco, &
                            zrefo, zrefdo, ztrao, ztrado, 1)



!$acc kernels loop    
do iplon = 1, ncol

! Combine clear and cloudy reflectivies and optical depths     

!$acc loop
      do iw = 1, 112
            
!$acc loop seq
            do jk=1,klev

! Combine clear and cloudy contributions for total sky
               !ikl = klev+1-jk 

! Direct beam transmittance        

               ze1 = (ztauo(iplon,iw,jk))  / prmu0(iplon)      
               zdbtmc = exp(-ze1)
               zdbt(iplon,iw,jk)   = zdbtmc
               ztdbt(iplon,iw,jk+1)   = zdbt(iplon,iw,jk)  *ztdbt(iplon,iw,jk)  

            enddo          
        end do
end do
!$acc end kernels

! compute the fluxes from the optical depths and reflectivities

                 
! Vertical quadrature for clear-sky fluxes

!$acc kernels 
do iplon = 1, ncol


! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

    do iw = 1, 112
          jb = ngb(iw)
          ibm = jb-15

            zgco(iplon,iw,klev+1)   =palbp(iplon,ibm) 
            zomco(iplon,iw,klev+1)  =palbd(iplon,ibm) 
    
    end do
end do
!$acc end kernels  

            call vrtqdr_sw(ncol, klev, &
                           zrefo  , zrefdo  , ztrao  , ztrado  , &
                           zdbt , zrdnd  , zgco, zomco, ztdbt  , &
                           zfd , zfu  , ztra)
            
            
            
           
            
! perform band integration for clear cases      
!$acc kernels loop
do iplon = 1, ncol
    
!$acc loop    
    do ikl=1,klev+1
            
      !$acc loop seq
      do iw = 1, 112
          jb = ngb(iw)
      
          jk=klev+2-ikl
          ibm = jb-15

! Apply adjustment for correct Earth/Sun distance and zenith angle to incoming solar flux
! No solar variability and no solar cycle
           if (isolvar .lt. 0) then
              zincflx = adjflux(jb)  * zsflxzen(iplon,iw)   * prmu0(iplon)           
           endif
! Solar variability with averaged or specified solar cycle
           if (isolvar .ge. 0) then
              zincflx = adjflux(jb)  * ssi(iplon,iw)   * prmu0(iplon)           
           endif


! Accumulate spectral fluxes over whole spectrum  
              
               pbbcu(iplon,ikl)  = pbbcu(iplon,ikl)  + zincflx*zfu(iplon,iw,jk)  
               pbbcd(iplon,ikl)  = pbbcd(iplon,ikl)  + zincflx*zfd(iplon,iw,jk)  
               pbbcddir(iplon,ikl)  = pbbcddir(iplon,ikl)  + zincflx*ztdbt(iplon,iw,jk)  
              

! Accumulate direct fluxes for UV/visible bands
               if (ibm >= 10 .and. ibm <= 13) then
                  puvcd(iplon,ikl)  = puvcd(iplon,ikl)  + zincflx*zfd(iplon,iw,jk)  
                  puvcddir(iplon,ikl)  = puvcddir(iplon,ikl)  + zincflx*ztdbt(iplon,iw,jk)  
                 
! Accumulate direct fluxes for near-IR bands
               else if (ibm == 14 .or. ibm <= 9) then  
                  pnicd(iplon,ikl)  = pnicd(iplon,ikl)  + zincflx*zfd(iplon,iw,jk)  
                  pnicddir(iplon,ikl)  = pnicddir(iplon,ikl)  + zincflx*ztdbt(iplon,iw,jk)  

               endif 

            enddo    

! End loop on jb, spectral band
      enddo

! End of longitude loop    
enddo               
!$acc end kernels


          
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  END CLEAR  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (cc==2) then

!$acc kernels 
do iplon = 1, ncol
    do iw = 1, 112
            do jk=1,klev

               ikl=klev+1-jk
               jb = ngb(iw)
               ibm = jb-15

               ze1 = ztaur(iplon,ikl,iw) + ptaua(iplon,ikl,ibm) * pomga(iplon, ikl, ibm) 
               ze2 = pasya(iplon, ikl, ibm) * pomga(iplon, ikl, ibm) * ptaua(iplon, ikl, ibm) / ze1
               ze1 = ze1/ (ztaur(iplon,ikl,iw)  + ztaug(iplon,ikl,iw) + ptaua(iplon,ikl,ibm)  )
               
               ! delta scale 
               zf = ze2*ze2
               zwf = ze1*zf
               ze1 = (ze1 - zwf) / (1.0 - zwf)
               ze2 = (ze2 - zf) / (1.0 - zf)
               
               
               ! delta scale
               zomco(iplon,iw,jk)   = (ztauo(iplon,iw,jk) * ze1  + ptaucmc(iplon,ikl,iw)  * pomgcmc(iplon,ikl,iw))
                     
               
               zgco(iplon, iw, jk) =  (ptaucmc(iplon,ikl,iw)  * pomgcmc(iplon,ikl,iw)  * pasycmc(iplon,ikl,iw) ) + &
                                      (ztauo(iplon, iw, jk) * ze1 * ze2)
               
               ztauo(iplon,iw,jk)   = ztauo(iplon,iw,jk) + ptaucmc(iplon,ikl,iw) 
     

               zgco(iplon,iw,jk)   = zgco(iplon, iw, jk) / zomco(iplon, iw, jk)
               zomco(iplon,iw,jk)  = zomco(iplon,iw,jk) / ztauo(iplon,iw,jk)
               
             
            enddo    
      end do
end do
!$acc end kernels



! Total sky reflectivities      
            call reftra_sw (ncol, nlayers, &
                            pcldfmc, zgco, prmu0, ztauo, zomco, &
                            zref, zrefd, ztra, ztrad, 0)
            



klev = nlayers




!$acc kernels loop    
do iplon = 1, ncol

!$acc loop
      do iw = 1, 112
            
!$acc loop seq
            do jk=1,klev

! Combine clear and cloudy contributions for total sky
               ikl = klev+1-jk 
               zclear = 1.0  - pcldfmc(iplon,ikl,iw) 
               zcloud = pcldfmc(iplon,ikl,iw) 

               zref(iplon,iw,jk)   = zclear*zrefo(iplon,iw,jk)   + zcloud*zref(iplon,iw,jk)  
               zrefd(iplon,iw,jk)  = zclear*zrefdo(iplon,iw,jk)   + zcloud*zrefd(iplon,iw,jk)  
               ztra(iplon,iw,jk)   = zclear*ztrao(iplon,iw,jk)   + zcloud*ztra(iplon,iw,jk)  
               ztrad(iplon,iw,jk)  = zclear*ztrado(iplon,iw,jk)   + zcloud*ztrad(iplon,iw,jk)  

! Clear + Cloud

               ze1 = ztauo(iplon,iw,jk )   / prmu0(iplon)   
               zdbtmo = exp(-ze1)            
               ze1 = (ztauo(iplon,iw,jk) - ptaucmc(iplon,ikl,iw))  / prmu0(iplon)           
               zdbtmc = exp(-ze1)

               zdbt(iplon,iw,jk)   = zclear*zdbtmc + zcloud*zdbtmo
               ztdbt(iplon,iw,jk+1)   = zdbt(iplon,iw,jk)  *ztdbt(iplon,iw,jk)  
            enddo          
        end do
end do
!$acc end kernels

!$acc kernels
zrdnd = 0.0
zgco = 0.0
zomco = 0.0
zfd = 0.0
zfu = 0.0
!$acc end kernels


!$acc kernels 
do iplon = 1, ncol

        
! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

    do iw = 1, 112
          jb = ngb(iw)
          ibm = jb-15

            zgco(iplon,iw,klev+1)   =palbp(iplon,ibm) 
            zomco(iplon,iw,klev+1)  =palbd(iplon,ibm) 
    
    end do
            enddo           
!$acc end kernels  
                 

      
! Vertical quadrature for cloudy fluxes


            call vrtqdr_sw(ncol, klev, &
                           zref, zrefd, ztra, ztrad, &
                           zdbt , zrdnd  , zgco, zomco , ztdbt  , &
                           zfd , zfu  ,  ztrao)

            

! Upwelling and downwelling fluxes at levels
!   Two-stream calculations go from top to bottom; 
!   layer indexing is reversed to go bottom to top for output arrays



  klev = nlayers
  

repclc = 1.e-12 

!$acc kernels loop
do iplon = 1, ncol
    
!$acc loop    
    do ikl=1,klev+1
            
      !$acc loop seq
      do iw = 1, 112
          jb = ngb(iw)
      
          jk=klev+2-ikl
          ibm = jb-15

! Apply adjustment for correct Earth/Sun distance and zenith angle to incoming solar flux
! No solar variability and no solar cycle
           if (isolvar .lt. 0) then
              zincflx = adjflux(jb)  * zsflxzen(iplon,iw)   * prmu0(iplon)           
           endif
! Solar variability with averaged or specified solar cycle
           if (isolvar .ge. 0) then
              zincflx = adjflux(jb)  * ssi(iplon,iw)   * prmu0(iplon)           
           endif



! Accumulate spectral fluxes over whole spectrum  
               pbbfu(iplon,ikl)  = pbbfu(iplon,ikl)  + zincflx*zfu(iplon,iw,jk)  
               pbbfd(iplon,ikl)  = pbbfd(iplon,ikl)  + zincflx*zfd(iplon,iw,jk)              
               pbbfddir(iplon,ikl)  = pbbfddir(iplon,ikl)  + zincflx*ztdbt(iplon,iw,jk)  

! Accumulate direct fluxes for UV/visible bands
               if (ibm >= 10 .and. ibm <= 13) then
                 
                  puvfd(iplon,ikl)  = puvfd(iplon,ikl)  + zincflx*zfd(iplon,iw,jk)  
                  puvfddir(iplon,ikl)  = puvfddir(iplon,ikl)  + zincflx*ztdbt(iplon,iw,jk)  
                 
                 
! Accumulate direct fluxes for near-IR bands
               else if (ibm == 14 .or. ibm <= 9) then  
                
                  pnifd(iplon,ikl)  = pnifd(iplon,ikl)  + zincflx*zfd(iplon,iw,jk)  
                  pnifddir(iplon,ikl)  = pnifddir(iplon,ikl)  + zincflx*ztdbt(iplon,iw,jk)  
                   
                 
               endif

            enddo

! End loop on jb, spectral band
         enddo             

! End of longitude loop    
enddo               
!$acc end kernels




else
!$acc kernels
    pbbfd = pbbcd
    pbbfu = pbbcu
    puvfd = puvcd
    puvfddir = puvcddir
    pnifd = pnicd
    pnifddir = pnicddir
!$acc end kernels    
end if


!$acc kernels
do iplon = 1, ncol
    
    do iw = 1, 112
        jb = ngb(iw)
        ibm = jb - 15
! Apply adjustment for correct Earth/Sun distance and zenith angle to incoming solar flux
! No solar variability and no solar cycle
           if (isolvar .lt. 0) then
              zincflx = adjflux(jb)  * zsflxzen(iplon,iw)   * prmu0(iplon)           
           endif
! Solar variability with averaged or specified solar cycle
           if (isolvar .ge. 0) then
              zincflx = adjflux(jb)  * ssi(iplon,iw)   * prmu0(iplon)           
           endif
            
         ! Accumulate surface direct fluxes for NIR
            if (ibm == 14 .or. ibm <= 8) then
               znirr(iplon) = znirr(iplon) + zincflx*ztdbt(iplon, iw, klev+1)        ! Direct flux
               znirf(iplon) = znirf(iplon) + zincflx*zfd(iplon, iw, klev+1)       ! Total flux
! Accumulate surface direct fluxes for PAR
            else if (ibm >= 10 .and. ibm <= 11) then
               zparr(iplon) = zparr(iplon) + zincflx*ztdbt(iplon, iw, klev+1)     ! Direct flux
               zparf(iplon) = zparf(iplon) + zincflx*zfd(iplon, iw, klev+1)       ! Total flux
! Accumulate surface direct fluxes for UV
            else if (ibm >= 12 .and. ibm <= 13) then
               zuvrr(iplon) = zuvrr(iplon) + zincflx*ztdbt(iplon, iw, klev+1)       ! Direct flux
               zuvrf(iplon) = zuvrf(iplon) + zincflx*zfd(iplon, iw, klev+1)      ! Total flux
            else if ( ibm==9) then
               zparr(iplon) = zparr(iplon) + 0.5*zincflx*ztdbt(iplon, iw, klev+1)     ! Direct flux
               zparf(iplon) = zparf(iplon) + 0.5*zincflx*zfd  (iplon, iw, klev+1)     ! Total flux
               znirr(iplon) = znirr(iplon) + 0.5*zincflx*ztdbt(iplon, iw, klev+1)     ! Direct flux
               znirf(iplon) = znirf(iplon) + 0.5*zincflx*zfd  (iplon, iw, klev+1)     ! Total flux
            endif
    end do
      enddo                    
!$acc end kernels



!!$acc end data



      end subroutine spcvmc_sw

             ! --------------------------------------------------------------------
      subroutine reftra_sw(ncol, nlayers, pcldfmc, pgg, prmuzl, ptau, pw, &
                           pref, prefd, ptra, ptrad, ac)
! --------------------------------------------------------------------
  
! Purpose: computes the reflectivity and transmissivity of a clear or 
!   cloudy layer using a choice of various approximations.
!
! Interface:  *rrtmg_sw_reftra* is called by *rrtmg_sw_spcvrt*
!
! Description:
! explicit arguments :
! --------------------
! inputs
! ------ 
!      lrtchk  = .t. for all layers in clear profile
!      lrtchk  = .t. for cloudy layers in cloud profile 
!              = .f. for clear layers in cloud profile
!      pgg     = assymetry factor
!      prmuz   = cosine solar zenith angle
!      ptau    = optical thickness
!      pw      = single scattering albedo
!
! outputs
! -------
!      pref    : collimated beam reflectivity
!      prefd   : diffuse beam reflectivity 
!      ptra    : collimated beam transmissivity
!      ptrad   : diffuse beam transmissivity
!
!
! Method:
! -------
!      standard delta-eddington, p.i.f.m., or d.o.m. layer calculations.
!      kmodts  = 1 eddington (joseph et al., 1976)
!              = 2 pifm (zdunkowski et al., 1980)
!              = 3 discrete ordinates (liou, 1973)
!
!
! Modifications:
! --------------
! Original: J-JMorcrette, ECMWF, Feb 2003
! Revised for F90 reformatting: MJIacono, AER, Jul 2006
! Revised to add exponential lookup table: MJIacono, AER, Aug 2007
!
! ------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------

      integer , intent(in) :: nlayers
      integer , intent(in) :: ncol

      real,  intent(in) :: pcldfmc(:,:,:)                           ! Logical flag for reflectivity and
                                                               ! and transmissivity calculation; 
                                                               !   Dimensions: (:)

      real , intent(in) gpu_device :: pgg(:,:,:)                        ! asymmetry parameter
                                                               !   Dimensions: (:)
      real , intent(in) gpu_device :: ptau(:,:,:)                       ! optical depth
                                                               !   Dimensions: (:)
      real , intent(in) gpu_device :: pw(:,:,:)                         ! single scattering albedo 
                                                               !   Dimensions: (:)
      real ,  intent(in) :: prmuzl(:)                       ! cosine of solar zenith angle
      integer, intent(in) :: ac

! ------- Output -------

      real , intent(out) gpu_device :: pref(:,:,:)                    ! direct beam reflectivity
                                                               !   Dimensions: (:+1)
      real , intent(out) gpu_device :: prefd(:,:,:)                   ! diffuse beam reflectivity
                                                               !   Dimensions: (:+1)
      real , intent(out) gpu_device :: ptra(:,:,:)                    ! direct beam transmissivity
                                                               !   Dimensions: (:+1)
      real , intent(out) gpu_device :: ptrad(:,:,:)                   ! diffuse beam transmissivity
                                                               !   Dimensions: (:+1)


! ------- Local -------

      integer  :: jk, jl, kmodts
      integer  :: itind, iplon, iw

      real  :: tblind
      real  :: za, za1, za2
      real  :: zbeta, zdend, zdenr, zdent
      real  :: ze1, ze2, zem1, zem2, zemm, zep1, zep2
      real  :: zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4, zgt
      real  :: zr1, zr2, zr3, zr4, zr5
      real  :: zrk, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
      real  :: zsr3, zt1, zt2, zt3, zt4, zt5, zto1
      real  :: zw, zwcrit, zwo, prmuz

      real , parameter :: eps = 1.e-08 

      ! MAT These are 8-byte versions of zw, zg, and zwo. This is done
      ! MAT to avoid a divide-by-zero in the zwo calculation below. More
      ! MAT information below.
      ! MAT NOTE: This is not an official fix, just a patch to allow work
      ! MAT       for now.

      real*8 :: zw8, zg8, zwo8

!     ------------------------------------------------------------------

! Initialize

    

      zsr3=sqrt(3. )
      zwcrit=0.9999995 
      kmodts=2
      
!$acc kernels loop
 do iplon=1,ncol
!$acc loop
    do iw=1,112
!$acc loop private(zgamma1, zgamma2, zgamma3, zgamma4)
      do jk=1, nlayers
           prmuz = prmuzl(iplon)
         if ((.not.(pcldfmc(iplon,nlayers+1-jk,iw))  > 1.e-12) .and. ac==0  ) then
            pref(iplon,iw,jk)   =0. 
            ptra(iplon,iw,jk)   =1. 
            prefd(iplon,iw,jk)  =0. 
            ptrad(iplon,iw,jk)  =1. 
         else
            zto1=ptau(iplon,iw,jk)  
            zw  =pw(iplon,iw,jk)  
            zg  =pgg(iplon,iw,jk)    

            ! MAT Move zw and zg into 8-byte reals to avoid
            ! MAT divide-by-zero in zwo calculation below

            zw8 = zw
            zg8 = zg

! General two-stream expressions

            zg3= 3.  * zg
           
               zgamma1= (8.  - zw * (5.  + zg3)) * 0.25 
               zgamma2=  3.  *(zw * (1.  - zg )) * 0.25 
               zgamma3= (2.  - zg3 * prmuz ) * 0.25 
       
            zgamma4= 1.  - zgamma3
    
! Recompute original s.s.a. to test for conservative solution

            ! MAT The issue with this is as follows. A column occurs for
            ! MAT which zw = 0.0249478 and zg = 0.503158. If you
            ! MAT calculate the denominator of zwo, you get a value of
            ! MAT -0.000000064995 which, apparently, is then flushed to
            ! MAT zero by the compiler. Then zwo = zw / 0. and bam.
            ! MAT To avoid this loss of precision, we calculate zwo
            ! MAT in 8-byte reals.
            ! MAT
            ! MAT Original code
            ! MAT zwo= zw / (1.  - (1.  - zw) * (zg / (1.  - zg))**2)
            ! MAT
            ! MAT New code
            zwo8= zw8 / (1.0d0  - (1.0d0  - zw8) * (zg8 / (1.0d0  - zg8))**2)

            ! MAT Put zwo8 into a 4-byte real
            zwo = zwo8

            ! END MAT EDITS
    
            if (zwo >= zwcrit) then
! Conservative scattering

               za  = zgamma1 * prmuz 
               za1 = za - zgamma3
               zgt = zgamma1 * zto1
        
! Homogeneous reflectance and transmittance,
! collimated beam

               ze1 = min ( zto1 / prmuz , 500. )


               ze2 = exp(-ze1)
               pref(iplon,iw,jk)   = (zgt - za1 * (1.  - ze2)) / (1.  + zgt)
               ptra(iplon,iw,jk)   = 1.  - pref(iplon,iw,jk)  

! isotropic incidence

               prefd(iplon,iw,jk)   = zgt / (1.  + zgt)
               ptrad(iplon,iw,jk)   = 1.  - prefd(iplon,iw,jk)          

! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
! table returns a transmittance of 1.0.
               if (ze2 .eq. 1.0 ) then 
                  pref(iplon,iw,jk)   = 0.0 
                  ptra(iplon,iw,jk)   = 1.0 
                  prefd(iplon,iw,jk)   = 0.0 
                  ptrad(iplon,iw,jk)   = 1.0 
               endif

            else
! Non-conservative scattering

               za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
               za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
               zrk = sqrt ( zgamma1**2 - zgamma2**2)
               zrp = zrk * prmuz               
               zrp1 = 1.  + zrp
               zrm1 = 1.  - zrp
               zrk2 = 2.  * zrk
               zrpp = 1.  - zrp*zrp
               zrkg = zrk + zgamma1
               zr1  = zrm1 * (za2 + zrk * zgamma3)
               zr2  = zrp1 * (za2 - zrk * zgamma3)
               zr3  = zrk2 * (zgamma3 - za2 * prmuz )
               zr4  = zrpp * zrkg
               zr5  = zrpp * (zrk - zgamma1)
               zt1  = zrp1 * (za1 + zrk * zgamma4)
               zt2  = zrm1 * (za1 - zrk * zgamma4)
               zt3  = zrk2 * (zgamma4 + za1 * prmuz )
               zt4  = zr4
               zt5  = zr5

! mji - reformulated code to avoid potential floating point exceptions
!               zbeta = - zr5 / zr4
               zbeta = (zgamma1 - zrk) / zrkg
!!
        
! Homogeneous reflectance and transmittance

               ze1 = min ( zrk * zto1, 5. )
               ze2 = min ( zto1 / prmuz , 5. )

           
! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               if (ze1 .le. od_lo) then 
                  zem1 = 1.  - ze1 + 0.5  * ze1 * ze1
                  zep1 = 1.  / zem1
               else
                  
                  zem1 = exp(-ze1)
                  zep1 = 1.  / zem1
               endif
               if (ze2 .le. od_lo) then 
                  zem2 = 1.  - ze2 + 0.5  * ze2 * ze2
                  zep2 = 1.  / zem2
               else
                  zem2 = exp(-ze2)
                  zep2 = 1.  / zem2
               endif



               zdenr = zr4*zep1 + zr5*zem1
               zdent = zt4*zep1 + zt5*zem1
               if (zdenr .ge. -eps .and. zdenr .le. eps) then
                  pref(iplon,iw,jk)   = eps
                  ptra(iplon,iw,jk)   = zem2
               else 
                  pref(iplon,iw,jk)   = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
                  ptra(iplon,iw,jk)   = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent
               endif


! diffuse beam

               zemm = zem1*zem1
               zdend = 1.  / ( (1.  - zbeta*zemm ) * zrkg)
               prefd(iplon,iw,jk)   =  zgamma2 * (1.  - zemm) * zdend
               ptrad(iplon,iw,jk)   =  zrk2*zem1*zdend

            endif

         endif         

      end do  
    end do
   end do
!$acc end kernels
end subroutine reftra_sw
                           

                           
                           
! --------------------------------------------------------------------------
      subroutine vrtqdr_sw(ncol, klev, &
                           pref, prefd, ptra, ptrad, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
                           pfd, pfu, ztdn)
! --------------------------------------------------------------------------
 
! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
! 
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
!
!-----------------------------------------------------------------------

! ------- Declarations -------

! Input

      integer , intent (in) :: klev                   ! number of model layers
      integer , intent (in) :: ncol
    

      real , intent(in) gpu_device :: pref(:,:,:)                      ! direct beam reflectivity
                                                              !   Dimensions: (:+1)
      real , intent(in) gpu_device :: prefd(:,:,:)                     ! diffuse beam reflectivity
                                                              !   Dimensions: (:+1)
      real , intent(in) gpu_device :: ptra(:,:,:)                      ! direct beam transmissivity
                                                              !   Dimensions: (:+1)
      real , intent(in) gpu_device :: ptrad(:,:,:)                     ! diffuse beam transmissivity
                                                              !   Dimensions: (:+1)

      real , intent(in) gpu_device :: pdbt(:,:,:)  
                                                              !   Dimensions: (:+1)
      real , intent(in) gpu_device :: ptdbt(:,:,:)  
                                                              !   Dimensions: (:+1)

      real , intent(inout) gpu_device :: prdnd(:,:,:)  
                                                              !   Dimensions: (:+1)
      real , intent(inout) gpu_device :: prup(:,:,:)  
                                                              !   Dimensions: (:+1)
      real , intent(inout) gpu_device  :: prupd(:,:,:)  
                                                              !   Dimensions: (:+1)
      real, intent(inout) gpu_device :: ztdn(:,:,:)
                                                              
! Output
      real , intent(out) gpu_device  :: pfd(:,:,:)                    ! downwelling flux (W/m2)
                                                              !   Dimensions: (:+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle
      real , intent(inout) gpu_device  :: pfu(:,:,:)                    ! upwelling flux (W/m2)
                                                              !   Dimensions: (:+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle
    
     
! Local

      integer  :: ikp, ikx, jk, iplon, iw

      real  :: zreflect, zreflectj
     
      

! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------
                   
! Link lowest layer with surface
! this kernel has a lot of dependencies



!$acc kernels loop
do iplon = 1, ncol

    !$acc loop private(zreflect)
    do iw = 1, 112
      
        
      zreflect = 1.  / (1.  - prefd(iplon,iw,klev+1)   * prefd(iplon,iw,klev)  )
      prup(iplon,iw,klev)   = pref(iplon,iw,klev)   + (ptrad(iplon,iw,klev)   * &
                 ((ptra(iplon,iw,klev)   - pdbt(iplon,iw,klev)  ) * prefd(iplon,iw,klev+1)   + &
                   pdbt(iplon,iw,klev)   * pref(iplon,iw,klev+1)  )) * zreflect
      prupd(iplon,iw,klev)   = prefd(iplon,iw,klev)   + ptrad(iplon,iw,klev)   * ptrad(iplon,iw,klev)   * &
                    prefd(iplon,iw,klev+1)   * zreflect

    end do
end do
!$acc end kernels
      
! Pass from bottom to top 
!$acc kernels loop
do iplon = 1, ncol
    !$acc loop    
    do iw = 1, 112

      !$acc loop seq 
      do jk = 1,klev-1
         ikp = klev+1-jk                       
         ikx = ikp-1
         zreflectj = 1.  / (1.  -prupd(iplon,iw,ikp)   * prefd(iplon,iw,ikx)  )
         prup(iplon,iw,ikx)   = pref(iplon,iw,ikx)   + (ptrad(iplon,iw,ikx)   * &
                   ((ptra(iplon,iw,ikx)   - pdbt(iplon,iw,ikx)  ) * prupd(iplon,iw,ikp)   + &
                     pdbt(iplon,iw,ikx)   * prup(iplon,iw,ikp)  )) * zreflectj
         prupd(iplon,iw,ikx)   = prefd(iplon,iw,ikx)   + ptrad(iplon,iw,ikx)   * ptrad(iplon,iw,ikx)   * &
                      prupd(iplon,iw,ikp)   * zreflectj
      enddo
    end do
end do
!$acc end kernels

!$acc kernels loop
do iplon = 1, ncol
       !$acc loop
        do iw = 1, 112

! Upper boundary conditions

      ztdn(iplon, iw, 1) = 1. 
      prdnd(iplon,iw,1)   = 0. 
      ztdn(iplon, iw, 2) = ptra(iplon,iw,1)  
      prdnd(iplon,iw,2)   = prefd(iplon,iw,1)  
       end do
end do
!$acc end kernels      
      
!$acc kernels loop
do iplon = 1, ncol
    !$acc loop
    do iw = 1, 112

! Pass from top to bottom
      !$acc loop seq
      do jk = 2,klev
         ikp = jk+1
         zreflect = 1.  / (1.  - prefd(iplon,iw,jk)   * prdnd(iplon,iw,jk)  )
         ztdn(iplon, iw, ikp) = ptdbt(iplon,iw,jk)   * ptra(iplon,iw,jk)   + &
                    (ptrad(iplon,iw,jk)   * ((ztdn(iplon, iw, jk) - ptdbt(iplon,iw,jk)  ) + &
                     ptdbt(iplon,iw,jk)   * pref(iplon,iw,jk)   * prdnd(iplon,iw,jk)  )) * zreflect
         prdnd(iplon,iw,ikp)   = prefd(iplon,iw,jk)   + ptrad(iplon,iw,jk)   * ptrad(iplon,iw,jk)   * &
                      prdnd(iplon,iw,jk)   * zreflect
      enddo
    end do
end do
!$acc end kernels
    
! Up and down-welling fluxes at levels

!$acc kernels loop
do iplon = 1, ncol
    !$acc loop
    do iw = 1, 112
      !$acc loop 
      do jk = 1,klev+1
         zreflect = 1.  / (1.  - prdnd(iplon,iw,jk)   * prupd(iplon,iw,jk)  )
         pfu(iplon,iw,jk)   = (ptdbt(iplon,iw,jk)   * prup(iplon,iw,jk)   + &
                      (ztdn(iplon, iw, jk) - ptdbt(iplon,iw,jk)  ) * prupd(iplon,iw,jk)  ) * zreflect
         pfd(iplon,iw,jk)   = ptdbt(iplon,iw,jk)   + (ztdn(iplon, iw, jk) - ptdbt(iplon,iw,jk)  + &
                      ptdbt(iplon,iw,jk)   * prup(iplon,iw,jk)   * prdnd(iplon,iw,jk)  ) * zreflect
      enddo
    end do
end do
!$acc end kernels
      
end subroutine vrtqdr_sw

      end module rrtmg_sw_spcvmc



