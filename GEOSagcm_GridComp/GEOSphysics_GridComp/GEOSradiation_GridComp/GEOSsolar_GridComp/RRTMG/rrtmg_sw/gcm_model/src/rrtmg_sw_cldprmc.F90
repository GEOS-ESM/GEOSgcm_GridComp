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

   use parrrsw, only : ngptsw, jpband, jpb1, jpb2
   use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                        extice2, ssaice2, asyice2, &
                        extice3, ssaice3, asyice3, fdlice3, &
                        extice4, ssaice4, asyice4, &
                        abari, bbari, cbari, dbari, ebari, fbari
   use rrsw_wvn, only : ngb, icxa

   implicit none

contains

   ! ----------------------------------------------------------------------
   subroutine cldprmc_sw(ncol, nlayers, iceflag, liqflag, &
                         cldfmc, ciwpmc, clwpmc, reicmc, relqmc, &
                         taormc, taucmc, ssacmc, asmcmc)
   ! ----------------------------------------------------------------------

      ! Compute the cloud optical properties for each cloudy layer
      !   and g-point interval for use by the McICA method.  
      ! Note: Only liqflag=1/iceflag=1,2,3 are available;
      ! (Hu & Stamnes, Ebert and Curry, Key, and Fu) are implemented. 

      ! ------- Input -------

      integer, intent(in) :: ncol            ! number of gridcols
      integer, intent(in) :: nlayers         ! number of layers
      integer, intent(in) :: iceflag         ! see definitions
      integer, intent(in) :: liqflag         ! see definitions

      real, intent(in) :: cldfmc(:,:,:)          ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptsw,nlayers)
      real, intent(in) :: ciwpmc(:,:,:)          ! cloud ice water path [mcica]
                                                      !    Dimensions: (ngptsw,nlayers)
      real, intent(in) :: clwpmc(:,:,:)          ! cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptsw,nlayers)
      real, intent(in) :: relqmc(:,:)           ! cloud liquid particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
      real, intent(in) :: reicmc(:,:)           ! cloud ice particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
                                                      ! specific definition of reicmc depends on setting of iceflag:
                                                      ! iceflag = 0: (inactive)
                                                      !              
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]

      ! ------- Output -------

      real, intent(out) :: taucmc(:,:,:)       ! cloud optical depth (delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)
      real, intent(out) :: ssacmc(:,:,:)       ! single scattering albedo (delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)
      real, intent(out) :: asmcmc(:,:,:)       ! asymmetry parameter (delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)
      real, intent(out) :: taormc(:,:,:)         ! cloud optical depth (non-delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)

      ! ------- Local -------

      integer :: ib, lay, istr, index, icx, ig, iplon

      real, parameter :: eps = 1.e-06      ! epsilon
      real, parameter :: cldmin = 1.e-20   ! minimum value for cloud quantities

      real :: cwp                          ! total cloud water path
      real :: radliq                       ! cloud liquid droplet radius (microns)
      real :: radice                       ! cloud ice effective size (microns)

      real :: factor
      real :: fint

      real :: tauiceorig, scatice, ssaice, tauice
      real :: tauliqorig, scatliq, ssaliq, tauliq

      real :: fdelta
      real :: extcoice, gice, ssacoice, forwice
      real :: extcoliq, gliq, ssacoliq, forwliq

      !$acc kernels loop present(cldfmc, ciwpmc, clwpmc, relqmc, reicmc, taucmc, ssacmc, asmcmc, taormc)
      do iplon = 1, ncol

         !$acc loop 
         do lay = 1, nlayers

            !$acc loop private(fdelta,extcoice,gice,ssacoice,forwice,extcoliq,gliq,ssacoliq,forwliq)
            do ig = 1, ngptsw 

!?pmn: ordering is reverse of efficient
!?pmn: what if cloud not present where are outputs zeroed ???

               cwp = ciwpmc(iplon,lay,ig) + clwpmc(iplon,lay,ig)  
               if (cldfmc(iplon,lay,ig) >= cldmin .and. cwp >= cldmin) then

                  ! ---------------------------------------------------------
                  ! Calculation of absorption coefficients due to ice clouds.
                  ! ---------------------------------------------------------

                  if (ciwpmc(iplon,lay,ig) == 0.) then

                     extcoice = 0.
                     ssacoice = 0.
                     gice     = 0.
                     forwice  = 0.

                  elseif (iceflag == 1) then

                     ! Note: This option uses Ebert and Curry approach for all particle sizes similar to
                     ! CAM3 implementation, though this is somewhat unjustified for large ice particles
                   
                     ib = ngb(ig)
                     ib = icxa(ib)

                     radice = reicmc(iplon,lay) 
                     extcoice = (abari(ib) + bbari(ib)/radice)
                     ssacoice = 1. - cbari(ib) - dbari(ib) * radice
                     gice = ebari(ib) + fbari(ib) * radice
                     ! Check to ensure upper limit of gice is within physical limits for large particles
                     if (gice >= 1.) gice = 1. - eps
                     forwice = gice*gice
! Check to ensure all calculated quantities are within physical limits.

                  elseif (iceflag == 2) then
                     
                     ! Ice particle effective radius is limited to 5.0 to 131.0 um.

                     radice = reicmc(iplon,lay) 
                     factor = (radice - 2.) / 3. 
                     index = int(factor)
                     if (index == 43) index = 42
                     fint = factor - float(index)
                     ib = ngb(ig)
                     extcoice = extice2(index,ib) + &
                        fint * (extice2(index+1,ib) - extice2(index,ib))
                     ssacoice = ssaice2(index,ib) + &
                        fint * (ssaice2(index+1,ib) - ssaice2(index,ib))
                     gice = asyice2(index,ib) + &
                        fint * (asyice2(index+1,ib) - asyice2(index,ib))
                     forwice = gice*gice
! Check to ensure all calculated quantities are within physical limits.

                  elseif (iceflag == 3) then
                    
                     ! Ice particle generalized effective size is limited to 5.0 to 140.0 um.

                     radice = reicmc(iplon,lay) 
                     factor = (radice - 2.) / 3. 
                     index = int(factor)
                     if (index == 46) index = 45
                     fint = factor - float(index)
                     ib = ngb(ig)
                     extcoice = extice3(index,ib) + &
                        fint * (extice3(index+1,ib) - extice3(index,ib))
                     ssacoice = ssaice3(index,ib) + &
                        fint * (ssaice3(index+1,ib) - ssaice3(index,ib))
                     gice = asyice3(index,ib) + &
                        fint * (asyice3(index+1,ib) - asyice3(index,ib))
                     fdelta = fdlice3(index,ib) + &
                        fint * (fdlice3(index+1,ib) - fdlice3(index,ib))
                  
                     forwice = fdelta + 0.5 / ssacoice
                     ! See Fu 1996 p. 2067 
                     if (forwice .gt. gice) forwice = gice
! Check to ensure all calculated quantities are within physical limits.  

                  elseif (iceflag == 4) then

                     ! ice particle effective diameter is limited to 1.0 to 200.0 um.

                     radice = reicmc(iplon,lay) 
                     factor = radice
                     index = int(factor)
                     fint = factor - float(index)
                     ib = ngb(ig)
                     extcoice = extice4(index,ib) + &
                        fint * (extice4(index+1,ib) - extice4(index,ib))
                     ssacoice = ssaice4(index,ib) + &
                        fint * (ssaice4(index+1,ib) - ssaice4(index,ib))
                     gice = asyice4(index,ib) + &
                        fint * (asyice4(index+1,ib) - asyice4(index,ib))
                     forwice = gice*gice
! Check to ensure all calculated quantities are within physical limits.
               
                  endif

                  ! -----------------------------------------------------------
                  ! Calculation of absorption coefficients due to water clouds.
                  ! -----------------------------------------------------------

                  if (clwpmc(iplon,lay,ig) == 0.) then

                     extcoliq = 0.
                     ssacoliq = 0.
                     gliq     = 0.
                     forwliq  = 0.

                  elseif (liqflag == 1) then

                     radliq = relqmc(iplon,lay) 
                      
                     index = int(radliq - 1.5)
                     if (index == 0) index = 1
                     if (index == 58) index = 57
                     fint = radliq - 1.5 - float(index)
                     ib = ngb(ig)
                     extcoliq = extliq1(index,ib) + &
                        fint * (extliq1(index+1,ib) - extliq1(index,ib))
                     ssacoliq = ssaliq1(index,ib) + &
                        fint * (ssaliq1(index+1,ib) - ssaliq1(index,ib))
                     if (fint < 0. .and. ssacoliq > 1.) &
                        ssacoliq = ssaliq1(index,ib)
                     gliq = asyliq1(index,ib) + &
                        fint * (asyliq1(index+1,ib) - asyliq1(index,ib))
                     forwliq = gliq*gliq
! Check to ensure all calculated quantities are within physical limits.
                  
                  endif
   
                  tauliqorig = clwpmc(iplon,lay,ig) * extcoliq
                  tauiceorig = ciwpmc(iplon,lay,ig) * extcoice
                  taormc(iplon,lay,ig) = tauliqorig + tauiceorig

                  ssaliq = ssacoliq * (1. - forwliq) / (1. - forwliq * ssacoliq)
                  tauliq = (1. - forwliq * ssacoliq) * tauliqorig
                  ssaice = ssacoice * (1. - forwice) / (1. - forwice * ssacoice)
                  tauice = (1. - forwice * ssacoice) * tauiceorig

                  scatliq = ssaliq * tauliq
                  scatice = ssaice * tauice
                  taucmc(iplon,lay,ig) = tauliq + tauice

                  ! Ensure non-zero taucmc and scatice
                  if (taucmc(iplon,lay,ig) == 0.) taucmc(iplon,lay,ig) = cldmin
                  if (scatice == 0.) scatice = cldmin

                  ssacmc(iplon,lay,ig) = (scatliq + scatice) / taucmc(iplon,lay,ig)  
   
                  if (iceflag == 3) then
                     ! In accordance with the 1996 Fu paper, equation A.3, 
                     ! the moments for ice were calculated depending on whether using spheres
                     ! or hexagonal ice crystals.
                     ! Set asymetry parameter to first moment (istr=1)
                     istr = 1
                     asmcmc(iplon,lay,ig) = (1./(scatliq+scatice)) * &
                        (scatliq * (gliq**istr - forwliq) / (1. - forwliq) + &
                         scatice * ((gice-forwice) / (1. - forwice))**istr)

                  else 
                     ! This code is the standard method for delta-m scaling. 
                     ! Set asymetry parameter to first moment (istr=1)
                     istr = 1
                     asmcmc(iplon,lay,ig) = &
                          (scatliq * (gliq**istr - forwliq) / (1. - forwliq) + &
                           scatice * (gice**istr - forwice) / (1. - forwice))  &
                        / (scatliq + scatice)
                  endif 

               endif  ! cloud present

            enddo  ! g-points
         enddo  ! layers
      enddo  ! columns
      !$acc end kernels

   end subroutine cldprmc_sw

end module rrtmg_sw_cldprmc

