!?pmn Check to ensure all calculated quantities are within physical limits.
!?pmn use of cldmin needs improving

! space saving function-like macros for linear interpolation
! of a 2D-varaible in first and second arguments as follows ...
#define LIN2_ARG1(VAR,I,J,FINT) (VAR(I,J) + FINT * (VAR(I+1,J)-VAR(I,J)))
#define LIN2_ARG2(VAR,I,J,FINT) (VAR(I,J) + FINT * (VAR(I,J+1)-VAR(I,J)))

module rrtmg_sw_cldprmc

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

   use parrrsw, only : ngptsw
   use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                        extice2, ssaice2, asyice2, &
                        extice3, ssaice3, asyice3, fdlice3, &
                        extice4, ssaice4, asyice4, &
                        abari, bbari, cbari, dbari, ebari, fbari
   use rrsw_wvn, only : ngb, icxa

   implicit none

contains

   ! ----------------------------------------------------------------------
   subroutine cldprmc_sw(pncol, ncol, nlayers, iceflag, liqflag, &
                         cldymc, ciwpmc, clwpmc, reicmc, relqmc, &
                         taormc, taucmc, ssacmc, asmcmc)
   ! ----------------------------------------------------------------------

      ! Compute the cloud optical properties for each cloudy layer
      !   and g-point interval for use by the McICA method.  
      ! Note: Only liqflag=1/iceflag=1,2,3 are available;
      ! (Hu & Stamnes, Ebert and Curry, Key, and Fu) are implemented. 

      ! ------- Input -------

      integer, intent(in) :: pncol           ! Dimensioned num of gridcols
      integer, intent(in) :: ncol            ! Actual number of gridcols
      integer, intent(in) :: nlayers         ! number of layers
      integer, intent(in) :: iceflag         ! see definitions
      integer, intent(in) :: liqflag         ! see definitions

      logical, intent(in) :: cldymc (nlayers,ngptsw,pncol)  ! cloudy or not? [mcica]
      real,    intent(in) :: ciwpmc (nlayers,ngptsw,pncol)  ! cloud ice water path [mcica]
      real,    intent(in) :: clwpmc (nlayers,ngptsw,pncol)  ! cloud liq water path [mcica]
      real,    intent(in) :: relqmc (nlayers,pncol)         ! cloud liq ptle eff radius (um)
      real,    intent(in) :: reicmc (nlayers,pncol)         ! cloud ice ptle eff radius (um)

      ! Specific definition of reicmc depends on iceflag:
      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
      !                 r_ec range is limited to 13. to 130. um;
      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996),
      !                 r_k range is limited to 5. to 131. um;
      ! iceflag = 3: generalized effective size, dge, (Fu, 1996), 
      !                 dge range is limited to 5. to 140. um [dge = 1.0315 * r_ec].

      ! ------- Output -------

      ! Note: inout because ASSUMED that clear values already set by default externally.
      real, intent(inout) :: taucmc (nlayers,ngptsw,pncol)  ! cloud opt depth (delta scaled)
      real, intent(inout) :: ssacmc (nlayers,ngptsw,pncol)  ! single scat albedo (delta scaled)
      real, intent(inout) :: asmcmc (nlayers,ngptsw,pncol)  ! asymmetry param (delta scaled)
      real, intent(inout) :: taormc (nlayers,ngptsw,pncol)  ! cloud opt depth (non-delta scaled)

      ! ------- Local -------

      real, parameter :: eps = 1.e-06      ! small number cf 1.
      real, parameter :: cldmin = 1.e-20   ! minimum value for cloud quantities

      integer :: ib, lay, istr, ig, icol

      ! table lookup parameters
      integer :: index
      real :: factor, fint

      ! intermediates
      real :: radliq, radice
      real :: tauiceorig, scatice, ssaice, tauice
      real :: tauliqorig, scatliq, ssaliq, tauliq
      real :: fdelta

      logical :: cloudy (nlayers,ngptsw,pncol)

      real, dimension (nlayers,ngptsw,pncol) :: &
         extcoice, gice, ssacoice, forwice, &
         extcoliq, gliq, ssacoliq, forwliq

      ! --------------------------------
      ! Locations where do cloud optics.
      ! --------------------------------

      do icol = 1, ncol
         cloudy(:,:,icol) = cldymc(:,:,icol) .and. &
            (ciwpmc(:,:,icol) + clwpmc(:,:,icol)) >= cldmin
      end do
!?pmn: perhaps update cldymc here and make intent(inout)

      ! ---------------------------------------------------------
      ! Calculation of absorption coefficients due to ice clouds.
      ! ---------------------------------------------------------

      if (iceflag == 1) then

         ! Note: This option uses Ebert and Curry approach for all particle sizes
         ! similar to CAM3 implementation, though this is somewhat unjustified for
         ! large ice particles

         do icol = 1, ncol
            do ig = 1, ngptsw 
               ib = icxa(ngb(ig))
               do lay = 1, nlayers
                  if (cloudy(lay,ig,icol)) then

                     if (ciwpmc(lay,ig,icol) == 0.) then
                        extcoice(lay,ig,icol) = 0.
                        ssacoice(lay,ig,icol) = 0.
                        gice    (lay,ig,icol) = 0.
                        forwice (lay,ig,icol) = 0.
                     else
                        radice = reicmc(lay,icol) 
                        extcoice(lay,ig,icol) = abari(ib) + bbari(ib) / radice
                        ssacoice(lay,ig,icol) = 1. - cbari(ib) - dbari(ib) * radice
                        gice    (lay,ig,icol) = ebari(ib) + fbari(ib) * radice
                        ! Ensure gice within physical limits for large particles
                        if (gice(lay,ig,icol) >= 1.) gice(lay,ig,icol) = 1. - eps
!?pmn change to max (..., 1.-eps) when NZD
                        forwice (lay,ig,icol) = gice(lay,ig,icol) * gice(lay,ig,icol)
                     endif

                  endif  ! cloud present
               enddo  ! layers
            enddo  ! g-points
         enddo  ! columns
      
      elseif (iceflag == 2) then

         ! Ice particle effective radius is limited to 5.0 to 131.0 um.

         do icol = 1, ncol
            do ig = 1, ngptsw 
               ib = ngb(ig)
               do lay = 1, nlayers
                  if (cloudy(lay,ig,icol)) then

                     if (ciwpmc(lay,ig,icol) == 0.) then
                        extcoice(lay,ig,icol) = 0.
                        ssacoice(lay,ig,icol) = 0.
                        gice    (lay,ig,icol) = 0.
                        forwice (lay,ig,icol) = 0.
                     else
                        radice = reicmc(lay,icol) 
                        factor = (radice - 2.) / 3. 
                        index = int(factor)
                        if (index == 43) index = 42
                        fint = factor - float(index)
                        extcoice(lay,ig,icol) = LIN2_ARG1(extice2,index,ib,fint)
                        ssacoice(lay,ig,icol) = LIN2_ARG1(ssaice2,index,ib,fint)
                        gice    (lay,ig,icol) = LIN2_ARG1(asyice2,index,ib,fint)
                        forwice (lay,ig,icol) = gice(lay,ig,icol) * gice(lay,ig,icol)
                     endif

                  endif  ! cloud present
               enddo  ! layers
            enddo  ! g-points
         enddo  ! columns
      
      elseif (iceflag == 3) then
                    
         ! Ice particle generalized effective size is limited to 5. to 140. um.

         do icol = 1, ncol
            do ig = 1, ngptsw 
               ib = ngb(ig)
               do lay = 1, nlayers
                  if (cloudy(lay,ig,icol)) then

                     if (ciwpmc(lay,ig,icol) == 0.) then
                        extcoice(lay,ig,icol) = 0.
                        ssacoice(lay,ig,icol) = 0.
                        gice    (lay,ig,icol) = 0.
                        forwice (lay,ig,icol) = 0.
                     else
                        radice = reicmc(lay,icol) 
                        factor = (radice - 2.) / 3. 
                        index = int(factor)
                        if (index == 46) index = 45
                        fint = factor - float(index)
                        extcoice(lay,ig,icol) = LIN2_ARG1(extice3,index,ib,fint)
                        ssacoice(lay,ig,icol) = LIN2_ARG1(ssaice3,index,ib,fint)
                        gice    (lay,ig,icol) = LIN2_ARG1(asyice3,index,ib,fint)
                        fdelta                = LIN2_ARG1(fdlice3,index,ib,fint)
                        forwice (lay,ig,icol) = fdelta + 0.5 / ssacoice(lay,ig,icol)
                        ! See Fu 1996 p. 2067 
                        if (forwice(lay,ig,icol) > gice(lay,ig,icol)) &
                           forwice(lay,ig,icol) = gice(lay,ig,icol)
                     endif

                  endif  ! cloud present
               enddo  ! layers
            enddo  ! g-points
         enddo  ! columns

      elseif (iceflag == 4) then

         ! ice particle effective diameter is limited to 1. to 200. um.

         do icol = 1, ncol
            do ig = 1, ngptsw 
               ib = ngb(ig)
               do lay = 1, nlayers
                  if (cloudy(lay,ig,icol)) then

                     if (ciwpmc(lay,ig,icol) == 0.) then
                        extcoice(lay,ig,icol) = 0.
                        ssacoice(lay,ig,icol) = 0.
                        gice    (lay,ig,icol) = 0.
                        forwice (lay,ig,icol) = 0.
                     else
                        radice = reicmc(lay,icol) 
                        factor = radice
                        index = int(factor)
                        fint = factor - float(index)
                        extcoice(lay,ig,icol) = LIN2_ARG1(extice4,index,ib,fint)
                        ssacoice(lay,ig,icol) = LIN2_ARG1(ssaice4,index,ib,fint)
                        gice    (lay,ig,icol) = LIN2_ARG1(asyice4,index,ib,fint)
                        forwice (lay,ig,icol) = gice(lay,ig,icol) * gice(lay,ig,icol)
                     endif

                  endif  ! cloud present
               enddo  ! layers
            enddo  ! g-points
         enddo  ! columns

      endif  ! ice options

      ! -----------------------------------------------------------
      ! Calculation of absorption coefficients due to water clouds.
      ! -----------------------------------------------------------

      if (liqflag == 1) then

         do icol = 1, ncol
            do ig = 1, ngptsw 
               ib = ngb(ig)
               do lay = 1, nlayers
                  if (cloudy(lay,ig,icol)) then

                     if (clwpmc(lay,ig,icol) == 0.) then
                        extcoliq(lay,ig,icol) = 0.
                        ssacoliq(lay,ig,icol) = 0.
                        gliq    (lay,ig,icol) = 0.
                        forwliq (lay,ig,icol) = 0.
                     else
                        radliq = relqmc(lay,icol) 
                        index = int(radliq - 1.5)
                        if (index == 0) index = 1
                        if (index == 58) index = 57
                        fint = radliq - 1.5 - float(index)
                        extcoliq(lay,ig,icol) = LIN2_ARG1(extliq1,index,ib,fint)
                        ssacoliq(lay,ig,icol) = LIN2_ARG1(ssaliq1,index,ib,fint)
                        if (fint < 0. .and. ssacoliq(lay,ig,icol) > 1.) &
                           ssacoliq(lay,ig,icol) = ssaliq1(index,ib)
                        gliq    (lay,ig,icol) = LIN2_ARG1(asyliq1,index,ib,fint)
                        forwliq (lay,ig,icol) = gliq(lay,ig,icol) * gliq(lay,ig,icol)
                     endif

                  endif  ! cloud present
               enddo  ! layers
            enddo  ! g-points
         enddo  ! columns

      endif  ! liquid options

      ! --------------------------------------------------------------------
      ! Perform delta-scaling and Combine ice and liquid optical properties.
      ! --------------------------------------------------------------------

      do icol = 1, ncol
         do ig = 1, ngptsw 
            do lay = 1, nlayers
               if (cloudy(lay,ig,icol)) then

                  tauliqorig = clwpmc(lay,ig,icol) * extcoliq(lay,ig,icol)
                  tauiceorig = ciwpmc(lay,ig,icol) * extcoice(lay,ig,icol)
                  taormc(lay,ig,icol) = tauliqorig + tauiceorig

                  ssaliq = ssacoliq(lay,ig,icol) * (1. - forwliq(lay,ig,icol)) &
                           / (1. - forwliq(lay,ig,icol) * ssacoliq(lay,ig,icol))
                  ssaice = ssacoice(lay,ig,icol) * (1. - forwice(lay,ig,icol)) &
                           / (1. - forwice(lay,ig,icol) * ssacoice(lay,ig,icol))
                  tauliq = (1. - forwliq(lay,ig,icol) * ssacoliq(lay,ig,icol)) * tauliqorig
                  tauice = (1. - forwice(lay,ig,icol) * ssacoice(lay,ig,icol)) * tauiceorig

                  scatliq = ssaliq * tauliq
                  scatice = ssaice * tauice
                  taucmc(lay,ig,icol) = tauliq + tauice

                  ! Ensure non-zero taucmc and scatice
!? pmn because of normalization below?
!? pmn cldmin use below is klugy ... find a better way when can go NZD.
!? this usage gives ssacmc = cldmin/cldmin = 1 for no cloud
                  if (taucmc(lay,ig,icol) == 0.) taucmc(lay,ig,icol) = cldmin
                  if (scatice == 0.) scatice = cldmin

                  ssacmc(lay,ig,icol) = (scatliq + scatice) / taucmc(lay,ig,icol)  
   
                  if (iceflag == 3) then
                     ! In accordance with the 1996 Fu paper, equation A.3, 
                     ! the moments for ice were calculated depending on whether using spheres
                     ! or hexagonal ice crystals.
                     ! Set asymetry parameter to first moment (istr=1)
                     istr = 1
                     asmcmc(lay,ig,icol) = (1./(scatliq+scatice)) * &
                        (scatliq * (gliq(lay,ig,icol)**istr - forwliq(lay,ig,icol)) &
                                   / (1. - forwliq(lay,ig,icol)) + &
                         scatice * ((gice(lay,ig,icol)-forwice(lay,ig,icol)) &
                                   / (1. - forwice(lay,ig,icol)))**istr)

                  else 
                     ! This code is the standard method for delta-m scaling. 
                     ! Set asymetry parameter to first moment (istr=1)
                     istr = 1
                     asmcmc(lay,ig,icol) = &
                          (scatliq * (gliq(lay,ig,icol)**istr - forwliq(lay,ig,icol)) &
                                     / (1. - forwliq(lay,ig,icol)) + &
                           scatice * (gice(lay,ig,icol)**istr - forwice(lay,ig,icol)) &
                                     / (1. - forwice(lay,ig,icol)))  &
                          / (scatliq + scatice)
                  endif 

               endif  ! cloud present
            enddo  ! layers
         enddo  ! g-points
      enddo  ! columns

   end subroutine cldprmc_sw

end module rrtmg_sw_cldprmc

