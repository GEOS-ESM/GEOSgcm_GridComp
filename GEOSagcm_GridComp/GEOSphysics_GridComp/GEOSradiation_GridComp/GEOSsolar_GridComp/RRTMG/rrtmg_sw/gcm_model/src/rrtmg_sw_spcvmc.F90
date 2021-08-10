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
      palbd, palbp, &
      pcldfmc, ptaucmc, pasycmc, pomgcmc, ptaormc, &
      ptaua, pasya, pomga, prmu0, adjflux, &
      isolvar, svar_f, svar_s, svar_i, &
      svar_f_bnd, svar_s_bnd, svar_i_bnd, &
      laytrop, jp, jt, jt1, &
      colch4, colco2, colh2o, colmol, colo2, colo3, &
      fac00, fac01, fac10, fac11, &
      selffac, selffrac, indself, forfac, forfrac, indfor, &
      pbbfd, pbbfu, pbbcd, pbbcu, puvfd, puvcd, pnifd, pnicd, &
      pbbfddir, pbbcddir, puvfddir, puvcddir, pnifddir, pnicddir, &
      znirr, znirf, zparr, zparf, zuvrr, zuvrf)
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
      integer, intent(in) :: laytrop (tncol)

      integer, intent(in) :: jp  (nlayers,tncol) 
      integer, intent(in) :: jt  (nlayers,tncol) 
      integer, intent(in) :: jt1 (nlayers,tncol) 

      real, intent(in) :: adjflux(jpband)             ! Earth/Sun distance adjustment

      ! Solar variability
      integer, intent(in) :: isolvar                  ! Flag for solar variability method
      real, intent(in) :: svar_f                      ! Solar variability facular multiplier
      real, intent(in) :: svar_s                      ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                      ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd(jpband)          ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd(jpband)          ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd(jpband)          ! Solar variability baseline irradiance multiplier (by band)

      real, intent(in) :: palbd(nbndsw,tncol)         ! surface albedo (diffuse)
      real, intent(in) :: palbp(nbndsw,tncol)         ! surface albedo (direct)
      real, intent(in) :: prmu0(tncol)                ! cosine of solar zenith angle

      real, intent(in) :: pcldfmc(nlayers,ngptsw,tncol)  ! cloud fraction [mcica]
      real, intent(in) :: ptaucmc(nlayers,ngptsw,tncol)  ! cloud optical depth [mcica]
      real, intent(in) :: pasycmc(nlayers,ngptsw,tncol)  ! cloud asymmetry parameter [mcica]
      real, intent(in) :: pomgcmc(nlayers,ngptsw,tncol)  ! cloud single scattering albedo [mcica]
      real, intent(in) :: ptaormc(nlayers,ngptsw,tncol)  ! cloud optical depth, non-delta scaled [mcica]
   
      real, intent(in) :: ptaua(nlayers+1,nbndsw,tncol)  ! aerosol optical depth
      real, intent(in) :: pasya(nlayers+1,nbndsw,tncol)  ! aerosol asymmetry parameter
      real, intent(in) :: pomga(nlayers+1,nbndsw,tncol)  ! aerosol single scattering albedo
                                                               
      real, intent(in) :: colh2o (nlayers,tncol) 
      real, intent(in) :: colco2 (nlayers,tncol) 
      real, intent(in) :: colch4 (nlayers,tncol) 
      real, intent(in) :: colo3  (nlayers,tncol) 
      real, intent(in) :: colo2  (nlayers,tncol) 
      real, intent(in) :: colmol (nlayers,tncol) 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlayers,tncol)
      integer, intent(in) :: indfor   (nlayers,tncol)
      real,    intent(in) :: selffac  (nlayers,tncol)
      real,    intent(in) :: selffrac (nlayers,tncol)
      real,    intent(in) :: forfac   (nlayers,tncol)
      real,    intent(in) :: forfrac  (nlayers,tncol)

      ! pressure and temperature interpolation coefficients
      real,    intent(in),  dimension (nlayers,tncol) &
         :: fac00, fac01, fac10, fac11

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

      real :: zgco   (nlayers,ngptsw,tncol)
      real :: zomco  (nlayers,ngptsw,tncol)  
      real :: ztauo  (nlayers,ngptsw,tncol)  

      real :: zdbt   (nlayers,  ngptsw,tncol)
      real :: ztdbt  (nlayers+1,ngptsw,tncol)   

      real :: zfd    (nlayers+1,ngptsw,tncol)
      real :: zfu    (nlayers+1,ngptsw,tncol)   

      real :: zref   (nlayers+1,ngptsw,tncol)
      real :: zrefo  (nlayers+1,ngptsw,tncol)  
      real :: zrefd  (nlayers+1,ngptsw,tncol)
      real :: zrefdo (nlayers+1,ngptsw,tncol)  
      real :: ztra   (nlayers+1,ngptsw,tncol)
      real :: ztrao  (nlayers+1,ngptsw,tncol)  
      real :: ztrad  (nlayers+1,ngptsw,tncol)
      real :: ztrado (nlayers+1,ngptsw,tncol)  

      real :: ztaur  (tncol,nlayers,ngptsw)
      real :: ztaug  (tncol,nlayers,ngptsw) 

      real :: zsflxzen (tncol,ngptsw)
      real :: ssi      (tncol,ngptsw)

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
!     zsflxzen = 0.
!     ssi      = 0.
      znirr    = 0.
      znirf    = 0.
      zparr    = 0.
      zparf    = 0.
      zuvrr    = 0.
      zuvrf    = 0.
      klev     = nlayers

      ! Calculate the optical depths for gaseous absorption and Rayleigh scattering     
      call taumol_sw( &
         tncol, ncol, nlayers, &
         colh2o, colco2, colch4, colo2, colo3, colmol, &
         laytrop, jp, jt, jt1, fac00, fac01, fac10, fac11, &
         selffac, selffrac, indself, forfac, forfrac, indfor, &
         isolvar, svar_f, svar_s, svar_i, svar_f_bnd, svar_s_bnd, svar_i_bnd, &
         ssi, zsflxzen, ztaug, ztaur)

      ! Set fixed boundary values.
      ! The surface (jk=klev+1) zref* and ztra* never change from these values.
      ! The TOA ztdbt (jk=1) likewise  never changes.
      do icol = 1,ncol
         do iw = 1,ngptsw
            jb = ngb(iw) ! SW band: jb = 16:29
            ibm = jb-15  !   => ibm = 1:14

            ! TOA direct beam    
            ztdbt (1,iw,icol) = 1. 
    
            ! Clear-sky Surface values
            ztrao (klev+1,iw,icol) = 0. 
            ztrado(klev+1,iw,icol) = 0. 
            zrefo (klev+1,iw,icol) = palbp(ibm,icol) 
            zrefdo(klev+1,iw,icol) = palbd(ibm,icol) 
           
            ! Total sky Surface values
            ztra  (klev+1,iw,icol) = 0. 
            ztrad (klev+1,iw,icol) = 0. 
            zref  (klev+1,iw,icol) = palbp(ibm,icol) 
            zrefd (klev+1,iw,icol) = palbd(ibm,icol) 
         end do
      end do

      do icol = 1,ncol
         do iw = 1,ngptsw
            do jk = 1,klev

               ikl = klev+1-jk
               jb = ngb(iw)
               ibm = jb-15

               ! Clear-sky optical parameters including aerosols
               ztauo(jk,iw,icol) = ztaur(icol,ikl,iw) + ztaug(icol,ikl,iw) + ptaua(ikl,ibm,icol)      
               zomco(jk,iw,icol) = ztaur(icol,ikl,iw) + ptaua(ikl,ibm,icol) * pomga(ikl,ibm,icol)
               zgco (jk,iw,icol) = pasya(ikl,ibm,icol) * pomga(ikl,ibm,icol) * ptaua(ikl,ibm,icol) / zomco(jk,iw,icol)   
               zomco(jk,iw,icol) = zomco(jk,iw,icol) / ztauo(jk,iw,icol)
               
               zf = zgco(jk,iw,icol)
               zf = zf * zf
               zwf = zomco(jk,iw,icol) * zf
               ztauo(jk,iw,icol) = (1. - zwf) * ztauo(jk,iw,icol)
               zomco(jk,iw,icol) = (zomco(jk,iw,icol) - zwf) / (1. - zwf)
               zgco (jk,iw,icol) = (zgco (jk,iw,icol) - zf ) / (1. - zf )

            enddo    
         end do
      end do

      ! Clear sky reflectivities
      call reftra_sw (tncol, ncol, nlayers, &
                      pcldfmc, zgco, prmu0, ztauo, zomco, &
                      zrefo, zrefdo, ztrao, ztrado, 1)

!x pmn   ! Combine clear and cloudy reflectivies and optical depths     
!x pmn         ! Combine clear and cloudy contributions for total sky
      do icol = 1,ncol
         do iw = 1,ngptsw
            do jk = 1,klev

               ! Direct beam transmittance        
               ze1 = ztauo(jk,iw,icol) / prmu0(icol)      
               zdbtmc = exp(-ze1)
               zdbt(jk,iw,icol) = zdbtmc
               ztdbt(jk+1,iw,icol) = zdbt(jk,iw,icol) * ztdbt(jk,iw,icol)  

            enddo          
        end do
      end do

      ! compute the fluxes from the optical depths and reflectivities

      ! Vertical quadrature for clear-sky fluxes

      call vrtqdr_sw(tncol, ncol, klev, &
                     zrefo, zrefdo, ztrao, ztrado, &
                     zdbt, ztdbt, &
                     zfd, zfu)

      ! perform band integration for clear cases      
      do icol = 1,ncol
         do ikl=1,klev+1
            do iw = 1,ngptsw

               jk=klev+2-ikl
               jb = ngb(iw)
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
              
               pbbcu   (icol,ikl) = pbbcu   (icol,ikl) + zincflx * zfu  (jk,iw,icol)  
               pbbcd   (icol,ikl) = pbbcd   (icol,ikl) + zincflx * zfd  (jk,iw,icol)  
               pbbcddir(icol,ikl) = pbbcddir(icol,ikl) + zincflx * ztdbt(jk,iw,icol)  
              
               ! Accumulate direct fluxes for UV/visible bands
               if (ibm >= 10 .and. ibm <= 13) then
                  puvcd   (icol,ikl)  = puvcd   (icol,ikl)  + zincflx * zfd  (jk,iw,icol)  
                  puvcddir(icol,ikl)  = puvcddir(icol,ikl)  + zincflx * ztdbt(jk,iw,icol)  
                 
               ! Accumulate direct fluxes for near-IR bands
               else if (ibm == 14 .or. ibm <= 9) then  
                  pnicd   (icol,ikl) = pnicd   (icol,ikl) + zincflx * zfd  (jk,iw,icol)  
                  pnicddir(icol,ikl) = pnicddir(icol,ikl) + zincflx * ztdbt(jk,iw,icol)  

               endif 

            enddo  ! spectral loop
         enddo  ! layer loop
      enddo  ! column loop

      !!!!!!!!!!!!!!!!
      !! END CLEAR  !!
      !!!!!!!!!!!!!!!!

      if (cc == 2) then

         do icol = 1,ncol
            do iw = 1,ngptsw
               do jk = 1,klev

                  ikl=klev+1-jk
                  jb = ngb(iw)
                  ibm = jb-15

                  ze1 = ztaur(icol,ikl,iw) + ptaua(ikl,ibm,icol) * pomga(ikl,ibm,icol) 
                  ze2 = pasya(ikl,ibm,icol) * pomga(ikl,ibm,icol) * ptaua(ikl,ibm,icol) / ze1
                  ze1 = ze1 / (ztaur(icol,ikl,iw) + ztaug(icol,ikl,iw) + ptaua(ikl,ibm,icol))
               
                  ! delta scale 
                  zf = ze2*ze2
                  zwf = ze1*zf
                  ze1 = (ze1 - zwf) / (1. - zwf)
                  ze2 = (ze2 - zf) / (1. - zf)
               
                  ! delta scale
                  zomco(jk,iw,icol) = ztauo(jk,iw,icol) * ze1 + ptaucmc(ikl,iw,icol) * pomgcmc(ikl,iw,icol)
                        
                  zgco (jk,iw,icol) = ptaucmc(ikl,iw,icol) * pomgcmc(ikl,iw,icol) * pasycmc(ikl,iw,icol) + &
                                         ztauo(jk,iw,icol) * ze1 * ze2
               
                  ztauo(jk,iw,icol) = ztauo(jk,iw,icol) + ptaucmc(ikl,iw,icol) 
     
                  zgco (jk,iw,icol) = zgco (jk,iw,icol) / zomco(jk,iw,icol)
                  zomco(jk,iw,icol) = zomco(jk,iw,icol) / ztauo(jk,iw,icol)
               
               enddo    
            end do
         end do

         ! Total sky reflectivities      
         call reftra_sw (tncol, ncol, nlayers, &
                         pcldfmc, zgco, prmu0, ztauo, zomco, &
                         zref, zrefd, ztra, ztrad, 0)
!?pmn       
         klev = nlayers

         do icol = 1,ncol
            do iw = 1,ngptsw
               do jk = 1,klev
                  ikl = klev+1-jk 

                  ! Combine clear and cloudy contributions for total sky

                  zclear = 1. - pcldfmc(ikl,iw,icol) 
                  zcloud = pcldfmc(ikl,iw,icol) 

                  zref (jk,iw,icol) = zclear * zrefo (jk,iw,icol) + zcloud * zref (jk,iw,icol)  
                  zrefd(jk,iw,icol) = zclear * zrefdo(jk,iw,icol) + zcloud * zrefd(jk,iw,icol)  
                  ztra (jk,iw,icol) = zclear * ztrao (jk,iw,icol) + zcloud * ztra (jk,iw,icol)  
                  ztrad(jk,iw,icol) = zclear * ztrado(jk,iw,icol) + zcloud * ztrad(jk,iw,icol)  

                  ! Clear + Cloud

                  ze1 = ztauo(jk,iw,icol) / prmu0(icol)   
                  zdbtmo = exp(-ze1)            
                  ze1 = (ztauo(jk,iw,icol) - ptaucmc(ikl,iw,icol)) / prmu0(icol)           
                  zdbtmc = exp(-ze1)

                  zdbt(jk,iw,icol) = zclear * zdbtmc + zcloud * zdbtmo
                  ztdbt(jk+1,iw,icol) = zdbt(jk,iw,icol) * ztdbt(jk,iw,icol)  

               enddo          
            end do
         end do

         ! Vertical quadrature for cloudy fluxes

         call vrtqdr_sw(tncol, ncol, klev, &
                        zref, zrefd, ztra, ztrad, &
                        zdbt, ztdbt, &
                        zfd, zfu)

         ! Upwelling and downwelling fluxes at levels
         !   Two-stream calculations go from top to bottom; 
         !   layer indexing is reversed to go bottom to top for output arrays

!?pmn
         klev = nlayers

         do icol = 1,ncol
            do ikl=1,klev+1
               do iw = 1,ngptsw

                  jk=klev+2-ikl
                  jb = ngb(iw)
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
                  pbbfu   (icol,ikl)  = pbbfu   (icol,ikl) + zincflx * zfu  (jk,iw,icol)  
                  pbbfd   (icol,ikl)  = pbbfd   (icol,ikl) + zincflx * zfd  (jk,iw,icol)              
                  pbbfddir(icol,ikl)  = pbbfddir(icol,ikl) + zincflx * ztdbt(jk,iw,icol)  

                  ! Accumulate direct fluxes for UV/visible bands
                  if (ibm >= 10 .and. ibm <= 13) then
                 
                     puvfd   (icol,ikl) = puvfd   (icol,ikl) + zincflx * zfd  (jk,iw,icol)  
                     puvfddir(icol,ikl) = puvfddir(icol,ikl) + zincflx * ztdbt(jk,iw,icol)  
                 
                  ! Accumulate direct fluxes for near-IR bands
                  else if (ibm == 14 .or. ibm <= 9) then  
                
                     pnifd   (icol,ikl) = pnifd   (icol,ikl) + zincflx * zfd  (jk,iw,icol)  
                     pnifddir(icol,ikl) = pnifddir(icol,ikl) + zincflx * ztdbt(jk,iw,icol)  
                   
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

      do icol = 1,ncol
         do iw = 1,ngptsw
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
               znirr(icol) = znirr(icol) + zincflx * ztdbt(klev+1,iw,icol)  ! Direct flux
               znirf(icol) = znirf(icol) + zincflx * zfd  (klev+1,iw,icol)  ! Total flux
            ! Accumulate surface direct fluxes for PAR
            else if (ibm >= 10 .and. ibm <= 11) then
               zparr(icol) = zparr(icol) + zincflx * ztdbt(klev+1,iw,icol)  ! Direct flux
               zparf(icol) = zparf(icol) + zincflx * zfd  (klev+1,iw,icol)  ! Total flux
            ! Accumulate surface direct fluxes for UV
            else if (ibm >= 12 .and. ibm <= 13) then
               zuvrr(icol) = zuvrr(icol) + zincflx * ztdbt(klev+1,iw,icol)  ! Direct flux
               zuvrf(icol) = zuvrf(icol) + zincflx * zfd  (klev+1,iw,icol)  ! Total flux
            else if ( ibm==9) then
               zparr(icol) = zparr(icol) + 0.5 * zincflx * ztdbt(klev+1,iw,icol)  ! Direct flux
               zparf(icol) = zparf(icol) + 0.5 * zincflx * zfd  (klev+1,iw,icol)  ! Total flux
               znirr(icol) = znirr(icol) + 0.5 * zincflx * ztdbt(klev+1,iw,icol)  ! Direct flux
               znirf(icol) = znirf(icol) + 0.5 * zincflx * zfd  (klev+1,iw,icol)  ! Total flux
            endif

         end do
      enddo                    

   end subroutine spcvmc_sw

   ! --------------------------------------------------------------------
   subroutine reftra_sw(tncol, ncol, nlayers, &
                        pcldfmc, pgg, prmuzl, ptau, pw, &
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

      integer, intent (in) :: tncol                   ! dimensioned num of gridcols
      integer, intent (in) :: ncol                    ! actual number of gridcols
      integer, intent (in) :: nlayers


      real, intent(in) :: pcldfmc (nlayers,ngptsw,tncol)   ! cloud fraction
!?pmn redo this with logical mcica flag

      real, intent(in) :: pgg     (nlayers,ngptsw,tncol)   ! asymmetry parameter
      real, intent(in) :: ptau    (nlayers,ngptsw,tncol)   ! optical depth
      real, intent(in) :: pw      (nlayers,ngptsw,tncol)   ! single scattering albedo 
      real, intent(in) :: prmuzl                 (tncol)   ! cosine of solar zenith angle

      integer, intent(in) :: ac

      ! ------- Output -------

      real, intent(out) :: pref  (nlayers+1,ngptsw,tncol)  ! direct beam reflectivity
      real, intent(out) :: prefd (nlayers+1,ngptsw,tncol)  ! diffuse beam reflectivity
      real, intent(out) :: ptra  (nlayers+1,ngptsw,tncol)  ! direct beam transmissivity
      real, intent(out) :: ptrad (nlayers+1,ngptsw,tncol)  ! diffuse beam transmissivity
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

      zsr3 = sqrt(3.)
      zwcrit = 0.9999995 
      kmodts = 2
      
      do icol = 1,ncol
         do iw = 1,ngptsw
            do jk = 1,nlayers

               prmuz = prmuzl(icol)
               if (.not.(pcldfmc(nlayers+1-jk,iw,icol) > 1.e-12) .and. ac==0) then

                  pref (jk,iw,icol) = 0. 
                  ptra (jk,iw,icol) = 1. 
                  prefd(jk,iw,icol) = 0. 
                  ptrad(jk,iw,icol) = 1. 

               else

                  zto1 = ptau(jk,iw,icol)  
                  zw   = pw  (jk,iw,icol)  
                  zg   = pgg (jk,iw,icol)    

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
                     pref(jk,iw,icol) = (zgt - za1 * (1. - ze2)) / (1. + zgt)
                     ptra(jk,iw,icol) = 1. - pref(jk,iw,icol)  

                     ! isotropic incidence

                     prefd(jk,iw,icol) = zgt / (1. + zgt)
                     ptrad(jk,iw,icol) = 1. - prefd(jk,iw,icol)          

                     ! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
                     ! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
                     ! table returns a transmittance of 1.
                     if (ze2 == 1.) then 
                        pref (jk,iw,icol) = 0. 
                        ptra (jk,iw,icol) = 1. 
                        prefd(jk,iw,icol) = 0. 
                        ptrad(jk,iw,icol) = 1. 
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
                        pref(jk,iw,icol) = eps
                        ptra(jk,iw,icol) = zem2
                     else 
                        pref(jk,iw,icol) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
                        ptra(jk,iw,icol) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent
                     endif

                     ! diffuse beam

                     zemm = zem1*zem1
                     zdend = 1. / ((1. - zbeta * zemm) * zrkg)
                     prefd(jk,iw,icol) = zgamma2 * (1. - zemm) * zdend
                     ptrad(jk,iw,icol) = zrk2 * zem1 * zdend

                  endif
               endif         

            end do  
         end do
      end do

   end subroutine reftra_sw
                           

   ! --------------------------------------------------------------------------
   subroutine vrtqdr_sw(tncol, ncol, klev, &
                        pref, prefd, ptra, ptrad, &
                        pdbt, ptdbt, &
                        pfd, pfu)
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

      integer, intent (in) :: tncol                   ! dimensioned num of gridcols
      integer, intent (in) :: ncol                    ! actual number of gridcols
      integer, intent (in) :: klev                    ! number of model layers
    
      real, intent(in) :: pref (klev+1,ngptsw,tncol)  ! direct beam reflectivity
      real, intent(in) :: prefd(klev+1,ngptsw,tncol)  ! diffuse reflectivity
      real, intent(in) :: ptra (klev+1,ngptsw,tncol)  ! direct beam transmissivity
      real, intent(in) :: ptrad(klev+1,ngptsw,tncol)  ! diffuse transmissivity

      real, intent(in) :: pdbt (klev,  ngptsw,tncol)  ! lyr mean dir beam transmittance
      real, intent(in) :: ptdbt(klev+1,ngptsw,tncol)  ! total direct beam transmittance

      ! ----- Output -----
      ! unadjusted for earth/sun distance or zenith angle

      real, intent(out) :: pfd(klev+1,ngptsw,tncol)   ! downwelling flux (W/m2)
      real, intent(out) :: pfu(klev+1,ngptsw,tncol)   ! upwelling   flux (W/m2)
    
      ! ----- Local -----

      integer :: ikp, ikx, jk, icol, iw
      real :: zreflect, zreflectj

      real :: ztdn  (klev+1,ngptsw,tncol)
      real :: prup  (klev+1,ngptsw,tncol)  
      real :: prupd (klev+1,ngptsw,tncol)  
      real :: prdnd (klev+1,ngptsw,tncol)  
     
      do icol = 1,ncol
         do iw = 1,ngptsw
      
            ! The klev+1 level of prup/prupd require palbp/palbd which
            ! are already available in fixed klev+1 level of pref/prefd.
            prup (klev+1,iw,icol) = pref (klev+1,iw,icol)
            prupd(klev+1,iw,icol) = prefd(klev+1,iw,icol)

            ! Link lowest layer with surface
            zreflect = 1. / (1. - prefd(klev+1,iw,icol) * prefd(klev,iw,icol))
            prup(klev,iw,icol) = pref(klev,iw,icol) + (ptrad(klev,iw,icol) * &
               ((ptra(klev,iw,icol) - pdbt(klev,iw,icol)) * prefd(klev+1,iw,icol) + &
               pdbt(klev,iw,icol) * pref(klev+1,iw,icol))) * zreflect
            prupd(klev,iw,icol) = prefd(klev,iw,icol) + ptrad(klev,iw,icol) * ptrad(klev,iw,icol) * &
               prefd(klev+1,iw,icol) * zreflect

          end do
      end do
      
      ! Pass from bottom to top 
      do icol = 1,ncol
         do iw = 1,ngptsw
            do jk = 1,klev-1

               ikp = klev+1-jk                       
               ikx = ikp-1
               zreflectj = 1. / (1. - prupd(ikp,iw,icol) * prefd(ikx,iw,icol))
               prup(ikx,iw,icol) = pref(ikx,iw,icol) + (ptrad(ikx,iw,icol) * &
                  ((ptra(ikx,iw,icol) - pdbt(ikx,iw,icol)) * prupd(ikp,iw,icol) + &
                  pdbt(ikx,iw,icol) * prup(ikp,iw,icol))) * zreflectj
               prupd(ikx,iw,icol) = prefd(ikx,iw,icol) + ptrad(ikx,iw,icol) * ptrad(ikx,iw,icol) * &
                  prupd(ikp,iw,icol) * zreflectj
            enddo
         end do
      end do

      do icol = 1,ncol
         do iw = 1,ngptsw

            ! Upper boundary conditions

            ztdn (1,iw,icol) = 1. 
            prdnd(1,iw,icol) = 0. 
            ztdn (2,iw,icol) = ptra (1,iw,icol)  
            prdnd(2,iw,icol) = prefd(1,iw,icol)  
         end do
      end do
      
      do icol = 1,ncol
         do iw = 1,ngptsw

            ! Pass from top to bottom

            do jk = 2,klev
               ikp = jk+1

               zreflect = 1. / (1. - prefd(jk,iw,icol) * prdnd(jk,iw,icol))
               ztdn(ikp,iw,icol) = ptdbt(jk,iw,icol) * ptra(jk,iw,icol) + &
                  (ptrad(jk,iw,icol) * ((ztdn(jk,iw,icol) - ptdbt(jk,iw,icol)) + &
                  ptdbt(jk,iw,icol) * pref(jk,iw,icol) * prdnd(jk,iw,icol))) * zreflect
               prdnd(ikp,iw,icol) = prefd(jk,iw,icol) + ptrad(jk,iw,icol) * ptrad(jk,iw,icol) * &
                      prdnd(jk,iw,icol) * zreflect

            enddo
         end do
      end do
    
      ! Up and down-welling fluxes at levels

      do icol = 1,ncol
         do iw = 1,ngptsw
            do jk = 1,klev+1
               zreflect = 1. / (1. - prdnd(jk,iw,icol) * prupd(jk,iw,icol))
               pfu(jk,iw,icol) = (ptdbt(jk,iw,icol) * prup(jk,iw,icol) + &
                  (ztdn(jk,iw,icol) - ptdbt(jk,iw,icol)) * prupd(jk,iw,icol)) * zreflect
               pfd(jk,iw,icol) = ptdbt(jk,iw,icol) + (ztdn(jk,iw,icol) - ptdbt(jk,iw,icol) + &
                  ptdbt(jk,iw,icol) * prup(jk,iw,icol) * prdnd(jk,iw,icol)) * zreflect
            enddo
         end do
      end do
      
   end subroutine vrtqdr_sw

end module rrtmg_sw_spcvmc
