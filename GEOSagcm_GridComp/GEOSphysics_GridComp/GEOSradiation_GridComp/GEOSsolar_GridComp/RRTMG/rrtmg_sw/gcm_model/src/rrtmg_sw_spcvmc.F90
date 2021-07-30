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

   use parrrsw, only : nbndsw, ngptsw, mxmol, jpband
   use rrsw_tbl, only : tblint, bpade, od_lo, exp_tbl
   use rrsw_wvn, only : ngc, ngs, ngb
   use rrtmg_sw_taumol, only: taumol_sw

   implicit none

contains

   ! ---------------------------------------------------------------------------
   subroutine spcvmc_sw ( &
      cc, tncol, ncol, nlayers, istart, iend, &
      pavel, tavel, palbd, palbp, &
      pcldfmc, ptaucmc, pasycmc, pomgcmc, ptaormc, &
      ptaua, pasya, pomga, prmu0, coldry,  adjflux, &
      isolvar, svar_f, svar_s, svar_i, &
      svar_f_bnd, svar_s_bnd, svar_i_bnd, &
      laytrop, laylow, jp, jt, jt1, &
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

      ! ------- Input -------

      integer, intent(in) :: tncol, ncol, cc
      integer, intent(in) :: nlayers
      integer, intent(in) :: istart
      integer, intent(in) :: iend
      integer, intent(in) :: laytrop(:)
      integer, intent(in) :: laylow(:)

      integer, intent(in) :: indfor(:,:) 
                                                               !   Dimensions: (nlayers)
      integer, intent(in) :: indself(:,:) 
                                                               !   Dimensions: (nlayers)
      integer, intent(in) :: jp(:,:) 
                                                               !   Dimensions: (nlayers)
      integer, intent(in) :: jt(:,:) 
                                                               !   Dimensions: (nlayers)
      integer, intent(in) :: jt1(:,:) 
                                                               !   Dimensions: (nlayers)

      real, intent(in) :: pavel(:,:)                     ! layer pressure (hPa, mb) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: tavel(:,:)                     ! layer temperature (K)
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: coldry(:,:)                    ! dry air column density (mol/cm2)
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: colmol(:,:) 
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

      real, intent(in) :: palbd(:,:)                     ! surface albedo (diffuse)
                                                               !   Dimensions: (nbndsw)
      real, intent(in) :: palbp(:,:)                     ! surface albedo (direct)
                                                               !   Dimensions: (nbndsw)
      real, intent(in) :: prmu0(:)                       ! cosine of solar zenith angle
      real, intent(in) :: pcldfmc(:,:,:)                 ! cloud fraction [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real, intent(in) :: ptaucmc(:,:,:)                 ! cloud optical depth [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real, intent(in) :: pasycmc(:,:,:)                 ! cloud asymmetry parameter [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real, intent(in) :: pomgcmc(:,:,:)                 ! cloud single scattering albedo [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
      real, intent(in) :: ptaormc(:,:,:)                 ! cloud optical depth, non-delta scaled [mcica]
                                                               !   Dimensions: (nlayers,ngptsw)
   
                                                               !   Dimensions: (nlayers,ngptsw)
      real, intent(in) :: ptaua(:,:,:)                  ! aerosol optical depth
                                                               !   Dimensions: (nlayers,nbndsw)
      real, intent(in) :: pasya(:,:,:)                  ! aerosol asymmetry parameter
                                                               !   Dimensions: (nlayers,nbndsw)
      real, intent(in) :: pomga(:,:,:)                  ! aerosol single scattering albedo
                                                               !   Dimensions: (nlayers,nbndsw)

                                                               
                                                               
      real, intent(in) :: colh2o(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: colco2(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: colch4(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: co2mult(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: colo3(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: colo2(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: coln2o(:,:) 
                                                               !   Dimensions: (nlayers)

      real, intent(in) :: forfac(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: forfrac(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: selffac(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: selffrac(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: fac00(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: fac01(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: fac10(:,:) 
                                                               !   Dimensions: (nlayers)
      real, intent(in) :: fac11(:,:) 
                                                               !   Dimensions: (nlayers)

! pmn why inout?
      real, intent(inout) :: zgco(tncol,ngptsw,nlayers+1), zomco(tncol,ngptsw,nlayers+1)  
      real, intent(inout) :: zrdnd(tncol,ngptsw,nlayers+1) 
      real, intent(inout) :: zref(tncol,ngptsw,nlayers+1)  , zrefo(tncol,ngptsw,nlayers+1)  
      real, intent(inout) :: zrefd(tncol,ngptsw,nlayers+1)  , zrefdo(tncol,ngptsw,nlayers+1)  
      real, intent(inout) :: ztauo(tncol,ngptsw,nlayers)  
      real, intent(inout) :: zdbt(tncol,ngptsw,nlayers+1)  ,ztdbt(tncol,ngptsw,nlayers+1)   
      real, intent(inout) :: ztra(tncol,ngptsw,nlayers+1)  , ztrao(tncol,ngptsw,nlayers+1)  
      real, intent(inout) :: ztrad(tncol,ngptsw,nlayers+1)  , ztrado(tncol,ngptsw,nlayers+1)  
      real, intent(inout) :: zfd(tncol,ngptsw,nlayers+1)  , zfu(tncol,ngptsw,nlayers+1)   
      real, intent(inout) :: ztaur(tncol,nlayers,ngptsw), ztaug(tncol,nlayers,ngptsw) 
      real, intent(inout) :: zsflxzen(tncol,ngptsw)
      real, intent(inout) :: ssi(tncol,ngptsw)
   
      ! ------- Output -------
                                                               !   All Dimensions: (nlayers+1)
      real, intent(out) :: pbbcd(:,:) 
      real, intent(out) :: pbbcu(:,:) 
      real, intent(out) :: pbbfd(:,:) 
      real, intent(out) :: pbbfu(:,:) 
      real, intent(out) :: pbbfddir(:,:) 
      real, intent(out) :: pbbcddir(:,:) 

      real, intent(out) :: puvcd(:,:) 
      real, intent(out) :: puvfd(:,:) 
      real, intent(out) :: puvcddir(:,:) 
      real, intent(out) :: puvfddir(:,:) 

      real, intent(out) :: pnicd(:,:) 
      real, intent(out) :: pnifd(:,:) 
      real, intent(out) :: pnicddir(:,:) 
      real, intent(out) :: pnifddir(:,:) 
      
      real, intent(out), dimension(:) :: znirr,znirf,zparr,zparf,zuvrr,zuvrf

      ! ------- Local -------

      integer :: klev
      integer :: ibm, ikl, ikp, ikx
      integer :: iw, jb, jg, jl, jk

      integer :: itind

      real :: tblind, ze1
      real :: zclear, zcloud

      real :: zincflx, ze2
     
      real :: zdbtmc, zdbtmo, zf, zgw, zreflect
      real :: zwf, tauorig

      integer :: icol

      ! ------------------------------------------------------------------

      pbbcd    = 0. 
      pbbcu    = 0. 
      pbbfd    = 0. 
      pbbfu    = 0. 
      pbbcddir = 0. 
      pbbfddir = 0. 
      puvcd    = 0. 
      puvfd    = 0. 
      puvcddir = 0. 
      puvfddir = 0. 
      pnicd    = 0. 
      pnifd    = 0. 
      pnicddir = 0. 
      pnifddir = 0.
      zsflxzen = 0.
      ssi      = 0.
      znirr    = 0.
      znirf    = 0.
      zparr    = 0.
      zparf    = 0.
      zuvrr    = 0.
      zuvrf    = 0.
      klev     = nlayers

      ! Calculate the optical depths for gaseous absorption and Rayleigh scattering     
      call taumol_sw( &
         ncol, nlayers, &
         colh2o, colco2, colch4, colo2, colo3, colmol, &
         laytrop, jp, jt, jt1, fac00, fac01, fac10, fac11, &
         selffac, selffrac, indself, forfac, forfrac, indfor, &
         isolvar, svar_f, svar_s, svar_i, svar_f_bnd, svar_s_bnd, svar_i_bnd, &
         ssi, zsflxzen, ztaug, ztaur)

      do icol = 1, ncol

         ! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

         do iw = 1, 112
            jb = ngb(iw)
            ibm = jb-15

!?pmn       ! Clear-sky    
!?pmn       !   TOA direct beam    
           
            ! Cloudy-sky    
            !   Surface values
            ztrao (icol,iw,klev+1) = 0. 
            ztrado(icol,iw,klev+1) = 0. 
            zrefo (icol,iw,klev+1) = palbp(icol,ibm) 
            zrefdo(icol,iw,klev+1) = palbd(icol,ibm) 
           
            ! Total sky    
            !   TOA direct beam    
            ztdbt (icol,iw,1)      = 1. 
    
            !   Surface values
            zdbt  (icol,iw,klev+1) = 0. 
            ztra  (icol,iw,klev+1) = 0. 
            ztrad (icol,iw,klev+1) = 0. 
            zref  (icol,iw,klev+1) = palbp(icol,ibm) 
            zrefd (icol,iw,klev+1) = palbd(icol,ibm) 

         end do
      end do

      do icol = 1, ncol

         do iw = 1,112

            do jk=1,klev

               ikl = klev+1-jk
               jb = ngb(iw)
               ibm = jb-15

               ! Clear-sky optical parameters including aerosols
               ztauo(icol,iw,jk) = ztaur(icol,ikl,iw) + ztaug(icol,ikl,iw) + ptaua(icol,ikl,ibm)      
               zomco(icol,iw,jk) = ztaur(icol,ikl,iw) + ptaua(icol,ikl,ibm) * pomga(icol,ikl,ibm)
               zgco (icol,iw,jk) = pasya(icol,ikl,ibm) * pomga(icol,ikl,ibm) * ptaua(icol,ikl,ibm) / zomco(icol,iw,jk)   
               zomco(icol,iw,jk) = zomco(icol,iw,jk) / ztauo(icol,iw,jk)
               
               zf = zgco(icol,iw,jk)
               zf = zf * zf
               zwf = zomco(icol,iw,jk) * zf
               ztauo(icol,iw,jk) = (1. - zwf) * ztauo(icol,iw,jk)
               zomco(icol,iw,jk) = (zomco(icol,iw,jk) - zwf) / (1. - zwf)
               zgco (icol,iw,jk) = (zgco(icol,iw,jk) - zf) / (1. - zf)

            enddo    
         end do

      end do

      ! Clear sky reflectivities
      call reftra_sw (ncol, nlayers, &
                      pcldfmc, zgco, prmu0, ztauo, zomco, &
                      zrefo, zrefdo, ztrao, ztrado, 1)

      do icol = 1, ncol

         ! Combine clear and cloudy reflectivies and optical depths     

         do iw = 1, 112
            
            do jk=1,klev

               ! Combine clear and cloudy contributions for total sky

               ! Direct beam transmittance        

               ze1 = ztauo(icol,iw,jk) / prmu0(icol)      
               zdbtmc = exp(-ze1)
               zdbt(icol,iw,jk) = zdbtmc
               ztdbt(icol,iw,jk+1) = zdbt(icol,iw,jk) * ztdbt(icol,iw,jk)  

            enddo          
        end do
      end do

      ! compute the fluxes from the optical depths and reflectivities

      ! Vertical quadrature for clear-sky fluxes

      do icol = 1, ncol

         ! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

         do iw = 1, 112
            jb = ngb(iw)
            ibm = jb-15

            zgco(icol,iw,klev+1)   =palbp(icol,ibm) 
            zomco(icol,iw,klev+1)  =palbd(icol,ibm) 
    
         end do
      end do

      call vrtqdr_sw(ncol, klev, &
                     zrefo, zrefdo, ztrao, ztrado, &
                     zdbt, zrdnd, zgco, zomco, ztdbt, &
                     zfd, zfu, ztra)

      ! perform band integration for clear cases      
      do icol = 1, ncol
    
         do ikl=1,klev+1
            
            do iw = 1, 112
               jb = ngb(iw)
      
!?pmn
               jk=klev+2-ikl
               ibm = jb-15

               ! Apply adjustment for correct Earth/Sun distance and zenith angle to incoming solar flux

               ! No solar variability and no solar cycle
               if (isolvar .lt. 0) then
                  zincflx = adjflux(jb) * zsflxzen(icol,iw) * prmu0(icol)           
               endif
               ! Solar variability with averaged or specified solar cycle
               if (isolvar .ge. 0) then
                  zincflx = adjflux(jb) * ssi(icol,iw)      * prmu0(icol)           
               endif

               ! Accumulate spectral fluxes over whole spectrum  
              
               pbbcu   (icol,ikl) = pbbcu   (icol,ikl) + zincflx * zfu  (icol,iw,jk)  
               pbbcd   (icol,ikl) = pbbcd   (icol,ikl) + zincflx * zfd  (icol,iw,jk)  
               pbbcddir(icol,ikl) = pbbcddir(icol,ikl) + zincflx * ztdbt(icol,iw,jk)  
              
               ! Accumulate direct fluxes for UV/visible bands
               if (ibm >= 10 .and. ibm <= 13) then
                  puvcd   (icol,ikl)  = puvcd   (icol,ikl)  + zincflx * zfd  (icol,iw,jk)  
                  puvcddir(icol,ikl)  = puvcddir(icol,ikl)  + zincflx * ztdbt(icol,iw,jk)  
                 
               ! Accumulate direct fluxes for near-IR bands
               else if (ibm == 14 .or. ibm <= 9) then  
                  pnicd   (icol,ikl) = pnicd   (icol,ikl) + zincflx * zfd  (icol,iw,jk)  
                  pnicddir(icol,ikl) = pnicddir(icol,ikl) + zincflx * ztdbt(icol,iw,jk)  

               endif 

            enddo  ! spectral loop
         enddo  ! layer loop
      enddo  ! column loop

      !!!!!!!!!!!!!!!!
      !! END CLEAR  !!
      !!!!!!!!!!!!!!!!

      if (cc == 2) then

         do icol = 1, ncol
            do iw = 1, 112
               do jk=1,klev

                  ikl=klev+1-jk
                  jb = ngb(iw)
                  ibm = jb-15

                  ze1 = ztaur(icol,ikl,iw) + ptaua(icol,ikl,ibm) * pomga(icol,ikl,ibm) 
                  ze2 = pasya(icol,ikl,ibm) * pomga(icol,ikl,ibm) * ptaua(icol,ikl,ibm) / ze1
                  ze1 = ze1 / (ztaur(icol,ikl,iw) + ztaug(icol,ikl,iw) + ptaua(icol,ikl,ibm))
               
                  ! delta scale 
                  zf = ze2*ze2
                  zwf = ze1*zf
                  ze1 = (ze1 - zwf) / (1. - zwf)
                  ze2 = (ze2 - zf) / (1. - zf)
               
                  ! delta scale
                  zomco(icol,iw,jk) = ztauo(icol,iw,jk) * ze1 + ptaucmc(icol,ikl,iw) * pomgcmc(icol,ikl,iw)
                        
                  zgco (icol,iw,jk) = ptaucmc(icol,ikl,iw) * pomgcmc(icol,ikl,iw) * pasycmc(icol,ikl,iw) + &
                                         ztauo(icol,iw,jk) * ze1 * ze2
               
                  ztauo(icol,iw,jk) = ztauo(icol,iw,jk) + ptaucmc(icol,ikl,iw) 
     
                  zgco (icol,iw,jk) = zgco (icol,iw,jk) / zomco(icol,iw,jk)
                  zomco(icol,iw,jk) = zomco(icol,iw,jk) / ztauo(icol,iw,jk)
               
               enddo    
            end do
         end do

         ! Total sky reflectivities      
         call reftra_sw (ncol, nlayers, &
                         pcldfmc, zgco, prmu0, ztauo, zomco, &
                         zref, zrefd, ztra, ztrad, 0)
!?pmn       
         klev = nlayers

         do icol = 1,ncol
            do iw = 1,112
               do jk=1,klev
                  ikl = klev+1-jk 

                  ! Combine clear and cloudy contributions for total sky

                  zclear = 1. - pcldfmc(icol,ikl,iw) 
                  zcloud = pcldfmc(icol,ikl,iw) 

                  zref (icol,iw,jk) = zclear * zrefo (icol,iw,jk) + zcloud * zref (icol,iw,jk)  
                  zrefd(icol,iw,jk) = zclear * zrefdo(icol,iw,jk) + zcloud * zrefd(icol,iw,jk)  
                  ztra (icol,iw,jk) = zclear * ztrao (icol,iw,jk) + zcloud * ztra (icol,iw,jk)  
                  ztrad(icol,iw,jk) = zclear * ztrado(icol,iw,jk) + zcloud * ztrad(icol,iw,jk)  

                  ! Clear + Cloud

                  ze1 = ztauo(icol,iw,jk) / prmu0(icol)   
                  zdbtmo = exp(-ze1)            
                  ze1 = (ztauo(icol,iw,jk) - ptaucmc(icol,ikl,iw)) / prmu0(icol)           
                  zdbtmc = exp(-ze1)

                  zdbt(icol,iw,jk) = zclear * zdbtmc + zcloud * zdbtmo
                  ztdbt(icol,iw,jk+1) = zdbt(icol,iw,jk) * ztdbt(icol,iw,jk)  

               enddo          
            end do
         end do

         zrdnd = 0.
         zgco  = 0.
         zomco = 0.
         zfd   = 0.
         zfu   = 0.

         do icol = 1, ncol
        
            ! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

            do iw = 1,112
               jb = ngb(iw)
               ibm = jb-15

               zgco (icol,iw,klev+1) = palbp(icol,ibm) 
               zomco(icol,iw,klev+1) = palbd(icol,ibm) 
    
            end do
         enddo           
                 
         ! Vertical quadrature for cloudy fluxes

         call vrtqdr_sw(ncol, klev, &
                        zref, zrefd, ztra, ztrad, &
                        zdbt, zrdnd, zgco, zomco, ztdbt, &
                        zfd, zfu, ztrao)

         ! Upwelling and downwelling fluxes at levels
         !   Two-stream calculations go from top to bottom; 
         !   layer indexing is reversed to go bottom to top for output arrays

!?pmn
         klev = nlayers

         do icol = 1, ncol
    
            do ikl=1,klev+1
            
               do iw = 1, 112
                  jb = ngb(iw)
      
!?pmn order
                  jk=klev+2-ikl
                  ibm = jb-15

                  ! Apply adjustment for correct Earth/Sun distance and zenith angle to incoming solar flux
                  ! No solar variability and no solar cycle
                  if (isolvar .lt. 0) then
                     zincflx = adjflux(jb) * zsflxzen(icol,iw) * prmu0(icol)           
                  endif
                  ! Solar variability with averaged or specified solar cycle
                  if (isolvar .ge. 0) then
                     zincflx = adjflux(jb) * ssi     (icol,iw) * prmu0(icol)           
                  endif

                  ! Accumulate spectral fluxes over whole spectrum  
                  pbbfu   (icol,ikl)  = pbbfu   (icol,ikl) + zincflx * zfu  (icol,iw,jk)  
                  pbbfd   (icol,ikl)  = pbbfd   (icol,ikl) + zincflx * zfd  (icol,iw,jk)              
                  pbbfddir(icol,ikl)  = pbbfddir(icol,ikl) + zincflx * ztdbt(icol,iw,jk)  

                  ! Accumulate direct fluxes for UV/visible bands
                  if (ibm >= 10 .and. ibm <= 13) then
                 
                     puvfd   (icol,ikl) = puvfd   (icol,ikl) + zincflx * zfd  (icol,iw,jk)  
                     puvfddir(icol,ikl) = puvfddir(icol,ikl) + zincflx * ztdbt(icol,iw,jk)  
                 
                  ! Accumulate direct fluxes for near-IR bands
                  else if (ibm == 14 .or. ibm <= 9) then  
                
                     pnifd   (icol,ikl) = pnifd   (icol,ikl) + zincflx * zfd  (icol,iw,jk)  
                     pnifddir(icol,ikl) = pnifddir(icol,ikl) + zincflx * ztdbt(icol,iw,jk)  
                   
                  endif

               enddo  ! spectral loop
            enddo  ! layer loop
         enddo  ! column loop

      else  ! cc /= 2

         pbbfd    = pbbcd
         pbbfu    = pbbcu
         puvfd    = puvcd
         puvfddir = puvcddir
         pnifd    = pnicd
         pnifddir = pnicddir

      end if

      do icol = 1, ncol
         do iw = 1, 112
            jb = ngb(iw)
            ibm = jb - 15

            ! Apply adjustment for correct Earth/Sun distance and zenith angle to incoming solar flux

            ! No solar variability and no solar cycle
            if (isolvar .lt. 0) then
               zincflx = adjflux(jb) * zsflxzen(icol,iw) * prmu0(icol)           
            endif
            ! Solar variability with averaged or specified solar cycle
            if (isolvar .ge. 0) then
               zincflx = adjflux(jb) * ssi     (icol,iw) * prmu0(icol)           
            endif
            
            ! Accumulate surface direct fluxes for NIR
            if (ibm == 14 .or. ibm <= 8) then
               znirr(icol) = znirr(icol) + zincflx * ztdbt(icol,iw,klev+1)  ! Direct flux
               znirf(icol) = znirf(icol) + zincflx * zfd  (icol,iw,klev+1)  ! Total flux
            ! Accumulate surface direct fluxes for PAR
            else if (ibm >= 10 .and. ibm <= 11) then
               zparr(icol) = zparr(icol) + zincflx * ztdbt(icol,iw,klev+1)  ! Direct flux
               zparf(icol) = zparf(icol) + zincflx * zfd  (icol,iw,klev+1)  ! Total flux
            ! Accumulate surface direct fluxes for UV
            else if (ibm >= 12 .and. ibm <= 13) then
               zuvrr(icol) = zuvrr(icol) + zincflx * ztdbt(icol,iw,klev+1)  ! Direct flux
               zuvrf(icol) = zuvrf(icol) + zincflx * zfd  (icol,iw,klev+1)  ! Total flux
            else if ( ibm==9) then
               zparr(icol) = zparr(icol) + 0.5 * zincflx * ztdbt(icol,iw,klev+1)  ! Direct flux
               zparf(icol) = zparf(icol) + 0.5 * zincflx * zfd  (icol,iw,klev+1)  ! Total flux
               znirr(icol) = znirr(icol) + 0.5 * zincflx * ztdbt(icol,iw,klev+1)  ! Direct flux
               znirf(icol) = znirf(icol) + 0.5 * zincflx * zfd  (icol,iw,klev+1)  ! Total flux
            endif

         end do
      enddo                    

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

      ! ------- Input -------

      integer, intent(in) :: nlayers
      integer, intent(in) :: ncol

      real, intent(in) :: pcldfmc(:,:,:)                           ! Logical flag for reflectivity and
                                                               ! and transmissivity calculation; 
                                                               !   Dimensions: (:)

      real, intent(in) :: pgg(:,:,:)                        ! asymmetry parameter
                                                               !   Dimensions: (:)
      real, intent(in) :: ptau(:,:,:)                       ! optical depth
                                                               !   Dimensions: (:)
      real, intent(in) :: pw(:,:,:)                         ! single scattering albedo 
                                                               !   Dimensions: (:)
      real, intent(in) :: prmuzl(:)                       ! cosine of solar zenith angle
      integer, intent(in) :: ac

      ! ------- Output -------

      real, intent(out) :: pref(:,:,:)                    ! direct beam reflectivity
                                                               !   Dimensions: (:+1)
      real, intent(out) :: prefd(:,:,:)                   ! diffuse beam reflectivity
                                                               !   Dimensions: (:+1)
      real, intent(out) :: ptra(:,:,:)                    ! direct beam transmissivity
                                                               !   Dimensions: (:+1)
      real, intent(out) :: ptrad(:,:,:)                   ! diffuse beam transmissivity
                                                               !   Dimensions: (:+1)
      ! ------- Local -------

      integer :: jk, jl, kmodts
      integer :: itind, icol, iw

      real :: tblind
      real :: za, za1, za2
      real :: zbeta, zdend, zdenr, zdent
      real :: ze1, ze2, zem1, zem2, zemm, zep1, zep2
      real :: zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4, zgt
      real :: zr1, zr2, zr3, zr4, zr5
      real :: zrk, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
      real :: zsr3, zt1, zt2, zt3, zt4, zt5, zto1
      real :: zw, zwcrit, zwo, prmuz

      real, parameter :: eps = 1.e-08 

      ! MAT These are 8-byte versions of zw, zg, and zwo. This is done
      ! MAT to avoid a divide-by-zero in the zwo calculation below. More
      ! MAT information below.
      ! MAT NOTE: This is not an official fix, just a patch to allow work
      ! MAT       for now.

      real*8 :: zw8, zg8, zwo8

      ! ------------------------------------------------------------------

      zsr3 = sqrt(3. )
      zwcrit = 0.9999995 
      kmodts = 2
      
      do icol=1,ncol
         do iw=1,112
            do jk=1,nlayers

               prmuz = prmuzl(icol)
               if ((.not.(pcldfmc(icol,nlayers+1-jk,iw)) > 1.e-12) .and. ac==0) then

                  pref (icol,iw,jk) = 0. 
                  ptra (icol,iw,jk) = 1. 
                  prefd(icol,iw,jk) = 0. 
                  ptrad(icol,iw,jk) = 1. 

               else

                  zto1 = ptau(icol,iw,jk)  
                  zw   = pw  (icol,iw,jk)  
                  zg   = pgg (icol,iw,jk)    

                  ! MAT Move zw and zg into 8-byte reals to avoid
                  ! MAT divide-by-zero in zwo calculation below

                  zw8 = zw
                  zg8 = zg

                  ! General two-stream expressions

                  zg3 = 3. * zg
           
                  zgamma1 = (8. - zw * (5. + zg3)) * 0.25 
                  zgamma2 = 3. * (zw * (1. - zg)) * 0.25 
                  zgamma3 = (2. - zg3 * prmuz) * 0.25 
                  zgamma4 = 1. - zgamma3
    
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
                  zwo8 = zw8 / (1.d0 - (1.d0 - zw8) * (zg8 / (1.d0 - zg8))**2)

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

                     ze1 = min(zto1/prmuz, 500.)
                     ze2 = exp(-ze1)
                     pref(icol,iw,jk) = (zgt - za1 * (1. - ze2)) / (1. + zgt)
                     ptra(icol,iw,jk) = 1. - pref(icol,iw,jk)  

                     ! isotropic incidence

                     prefd(icol,iw,jk) = zgt / (1. + zgt)
                     ptrad(icol,iw,jk) = 1. - prefd(icol,iw,jk)          

                     ! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
                     ! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
                     ! table returns a transmittance of 1.
                     if (ze2 == 1.) then 
                        pref (icol,iw,jk) = 0. 
                        ptra (icol,iw,jk) = 1. 
                        prefd(icol,iw,jk) = 0. 
                        ptrad(icol,iw,jk) = 1. 
                     endif

                  else
                     ! Non-conservative scattering

                     za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
                     za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
                     zrk = sqrt(zgamma1**2 - zgamma2**2)
                     zrp = zrk * prmuz               
                     zrp1 = 1. + zrp
                     zrm1 = 1. - zrp
                     zrk2 = 2. * zrk
                     zrpp = 1. - zrp*zrp
                     zrkg = zrk + zgamma1
                     zr1  = zrm1 * (za2 + zrk * zgamma3)
                     zr2  = zrp1 * (za2 - zrk * zgamma3)
                     zr3  = zrk2 * (zgamma3 - za2 * prmuz)
                     zr4  = zrpp * zrkg
                     zr5  = zrpp * (zrk - zgamma1)
                     zt1  = zrp1 * (za1 + zrk * zgamma4)
                     zt2  = zrm1 * (za1 - zrk * zgamma4)
                     zt3  = zrk2 * (zgamma4 + za1 * prmuz)
                     zt4  = zr4
                     zt5  = zr5

                     ! mji - reformulated code to avoid potential floating point exceptions
                     !               zbeta = - zr5 / zr4
                     zbeta = (zgamma1 - zrk) / zrkg
        
                     ! Homogeneous reflectance and transmittance

                     ze1 = min(zrk * zto1, 5.)
                     ze2 = min(zto1 / prmuz, 5.)

!?pmn no LUT!
                     ! Use exponential lookup table for transmittance, or expansion of 
                     ! exponential for low tau
                     if (ze1 <= od_lo) then 
                        zem1 = 1. - ze1 + 0.5 * ze1 * ze1
                        zep1 = 1. / zem1
                     else
                        zem1 = exp(-ze1)
                        zep1 = 1. / zem1
                     endif
                     if (ze2 <= od_lo) then 
                        zem2 = 1. - ze2 + 0.5 * ze2 * ze2
                        zep2 = 1. / zem2
                     else
                        zem2 = exp(-ze2)
                        zep2 = 1. / zem2
                     endif

                     zdenr = zr4*zep1 + zr5*zem1
                     zdent = zt4*zep1 + zt5*zem1
                     if (zdenr .ge. -eps .and. zdenr .le. eps) then
                        pref(icol,iw,jk) = eps
                        ptra(icol,iw,jk) = zem2
                     else 
                        pref(icol,iw,jk) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
                        ptra(icol,iw,jk) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent
                     endif

                     ! diffuse beam

                     zemm = zem1*zem1
                     zdend = 1. / ((1. - zbeta * zemm) * zrkg)
                     prefd(icol,iw,jk) = zgamma2 * (1. - zemm) * zdend
                     ptrad(icol,iw,jk) = zrk2 * zem1 * zdend

                  endif
               endif         

            end do  
         end do
      end do

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

      ! ----- Input -----

      integer, intent (in) :: ncol                   ! number of gridcols
      integer, intent (in) :: klev                   ! number of model layers
    
      real, intent(in) :: pref(:,:,:)                      ! direct beam reflectivity
                                                              !   Dimensions: (:+1)
      real, intent(in) :: prefd(:,:,:)                     ! diffuse beam reflectivity
                                                              !   Dimensions: (:+1)
      real, intent(in) :: ptra(:,:,:)                      ! direct beam transmissivity
                                                              !   Dimensions: (:+1)
      real, intent(in) :: ptrad(:,:,:)                     ! diffuse beam transmissivity
                                                              !   Dimensions: (:+1)

      real, intent(in) :: pdbt(:,:,:)  
                                                              !   Dimensions: (:+1)
      real, intent(in) :: ptdbt(:,:,:)  
                                                              !   Dimensions: (:+1)

      real, intent(inout) :: prdnd(:,:,:)  
                                                              !   Dimensions: (:+1)
      real, intent(inout) :: prup(:,:,:)  
                                                              !   Dimensions: (:+1)
      real, intent(inout) :: prupd(:,:,:)  
                                                              !   Dimensions: (:+1)
      real, intent(inout) :: ztdn(:,:,:)
                                                              
      ! ----- Output -----

      real, intent(out) :: pfd(:,:,:)                    ! downwelling flux (W/m2)
                                                              !   Dimensions: (:+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle
      real, intent(inout) :: pfu(:,:,:)                    ! upwelling flux (W/m2)
                                                              !   Dimensions: (:+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle
    
      ! ----- Local -----

      integer :: ikp, ikx, jk, icol, iw
      real :: zreflect, zreflectj
     
      ! ----- Definitions -----
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

      do icol = 1,ncol
         do iw = 1,112
      
            zreflect = 1. / (1. - prefd(icol,iw,klev+1) * prefd(icol,iw,klev))
            prup(icol,iw,klev) = pref(icol,iw,klev) + (ptrad(icol,iw,klev) * &
               ((ptra(icol,iw,klev) - pdbt(icol,iw,klev)) * prefd(icol,iw,klev+1) + &
               pdbt(icol,iw,klev) * pref(icol,iw,klev+1))) * zreflect
            prupd(icol,iw,klev) = prefd(icol,iw,klev) + ptrad(icol,iw,klev) * ptrad(icol,iw,klev) * &
               prefd(icol,iw,klev+1) * zreflect

          end do
      end do
      
      ! Pass from bottom to top 
      do icol = 1, ncol
         do iw = 1, 112
            do jk = 1,klev-1

               ikp = klev+1-jk                       
               ikx = ikp-1
               zreflectj = 1. / (1. - prupd(icol,iw,ikp) * prefd(icol,iw,ikx))
               prup(icol,iw,ikx) = pref(icol,iw,ikx) + (ptrad(icol,iw,ikx) * &
                  ((ptra(icol,iw,ikx) - pdbt(icol,iw,ikx)) * prupd(icol,iw,ikp) + &
                  pdbt(icol,iw,ikx) * prup(icol,iw,ikp))) * zreflectj
               prupd(icol,iw,ikx) = prefd(icol,iw,ikx) + ptrad(icol,iw,ikx) * ptrad(icol,iw,ikx) * &
                  prupd(icol,iw,ikp) * zreflectj
            enddo
         end do
      end do

      do icol = 1, ncol
         do iw = 1, 112

            ! Upper boundary conditions

            ztdn (icol,iw,1) = 1. 
            prdnd(icol,iw,1) = 0. 
            ztdn (icol,iw,2) = ptra (icol,iw,1)  
            prdnd(icol,iw,2) = prefd(icol,iw,1)  
         end do
      end do
      
      do icol = 1,ncol
         do iw = 1,112

            ! Pass from top to bottom

            do jk = 2,klev
               ikp = jk+1

               zreflect = 1. / (1. - prefd(icol,iw,jk) * prdnd(icol,iw,jk))
               ztdn(icol,iw,ikp) = ptdbt(icol,iw,jk) * ptra(icol,iw,jk) + &
                  (ptrad(icol,iw,jk) * ((ztdn(icol,iw,jk) - ptdbt(icol,iw,jk)) + &
                  ptdbt(icol,iw,jk) * pref(icol,iw,jk) * prdnd(icol,iw,jk))) * zreflect
               prdnd(icol,iw,ikp) = prefd(icol,iw,jk) + ptrad(icol,iw,jk) * ptrad(icol,iw,jk) * &
                      prdnd(icol,iw,jk) * zreflect

            enddo
         end do
      end do
    
      ! Up and down-welling fluxes at levels

      do icol = 1,ncol
         do iw = 1,112
            do jk = 1,klev+1
               zreflect = 1. / (1. - prdnd(icol,iw,jk) * prupd(icol,iw,jk))
               pfu(icol,iw,jk) = (ptdbt(icol,iw,jk) * prup(icol,iw,jk) + &
                  (ztdn(icol,iw,jk) - ptdbt(icol,iw,jk)) * prupd(icol,iw,jk)) * zreflect
               pfd(icol,iw,jk) = ptdbt(icol,iw,jk) + (ztdn(icol,iw,jk) - ptdbt(icol,iw,jk) + &
                  ptdbt(icol,iw,jk) * prup(icol,iw,jk) * prdnd(icol,iw,jk)) * zreflect
            enddo
         end do
      end do
      
   end subroutine vrtqdr_sw

end module rrtmg_sw_spcvmc
