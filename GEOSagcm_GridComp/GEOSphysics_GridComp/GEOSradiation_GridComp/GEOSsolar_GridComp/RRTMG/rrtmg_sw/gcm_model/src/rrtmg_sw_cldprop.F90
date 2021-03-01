!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

      module rrtmg_sw_cldprop

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
      use parrrsw, only : nbndsw, jpband, jpb1, jpb2
      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                           extice2, ssaice2, asyice2, &
                           extice3, ssaice3, asyice3, fdlice3, &
                           abari, bbari, cbari, dbari, ebari, fbari
      use rrsw_wvn, only : wavenum1, wavenum2
      use rrsw_vsn, only : hvrcld, hnamcld

      implicit none

      contains

! ----------------------------------------------------------------------------
      subroutine cldprop_sw(nlayers, inflag, iceflag, liqflag, cldfrac, &
                            tauc, ssac, asmc, fsfc, ciwp, clwp, rei, rel, &
                            taucldorig, taucloud, ssacloud, asmcloud)
! ----------------------------------------------------------------------------

! Purpose: Compute the cloud optical properties for each cloudy layer.
! Note: Only inflag = 0 and inflag=2/liqflag=1/iceflag=1,2,3 are available;
! (Hu & Stamnes, Ebert and Curry, Key, and Fu) are implemented.

! ------- Input -------

      integer , intent(in) :: nlayers         ! total number of layers
      integer , intent(in) :: inflag          ! see definitions
      integer , intent(in) :: iceflag         ! see definitions
      integer , intent(in) :: liqflag         ! see definitions

      real , intent(in) :: cldfrac(:)         ! cloud fraction
                                                      !    Dimensions: (nlayers)
      real , intent(in) :: ciwp(:)            ! cloud ice water path
                                                      !    Dimensions: (nlayers)
      real , intent(in) :: clwp(:)            ! cloud liquid water path
                                                      !    Dimensions: (nlayers)
      real , intent(in) :: rei(:)             ! cloud ice particle effective size (microns)
                                                      !    Dimensions: (nlayers)
                                                      ! specific definition of rei depends on setting of iceflag:
                                                      ! iceflag = 0: (inactive)
                                                      !              
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]
      real , intent(in) :: rel(:)             ! cloud liquid particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
      real , intent(in) :: tauc(:,:)          ! cloud optical depth
                                                      !    Dimensions: (nbndsw,nlayers)
      real , intent(in) :: ssac(:,:)          ! single scattering albedo
                                                      !    Dimensions: (nbndsw,nlayers)
      real , intent(in) :: asmc(:,:)          ! asymmetry parameter
                                                      !    Dimensions: (nbndsw,nlayers)
      real , intent(in) :: fsfc(:,:)          ! forward scattering fraction
                                                      !    Dimensions: (nbndsw,nlayers)

! ------- Output -------

      real , intent(out) :: taucloud(:,:)     ! cloud optical depth (delta scaled)
                                                      !    Dimensions: (nlayers,jpband)
      real , intent(out) :: taucldorig(:,:)   ! cloud optical depth (non-delta scaled)
                                                      !    Dimensions: (nlayers,jpband)
      real , intent(out) :: ssacloud(:,:)     ! single scattering albedo (delta scaled)
                                                      !    Dimensions: (nlayers,jpband)
      real , intent(out) :: asmcloud(:,:)     ! asymmetry parameter (delta scaled)
                                                      !    Dimensions: (nlayers,jpband)

! ------- Local -------

!      integer  :: ncbands
      integer  :: ib, ib1, ib2, lay, istr, index, icx

      real , parameter :: eps = 1.e-06      ! epsilon
      real , parameter :: cldmin = 1.e-20   ! minimum value for cloud quantities
      real  :: cwp                            ! total cloud water path
      real  :: radliq                         ! cloud liquid droplet radius (microns)
      real  :: radice                         ! cloud ice effective size (microns)
      real  :: factor
      real  :: fint
      real  :: tauctot(nlayers)               ! band integrated cloud optical depth

      real  :: taucldorig_a, ssacloud_a, taucloud_a, ffp, ffp1, ffpssa
      real  :: tauiceorig, scatice, ssaice, tauice, tauliqorig, scatliq, ssaliq, tauliq

      real  :: fdelta(jpb1:jpb2)
      real  :: extcoice(jpb1:jpb2), gice(jpb1:jpb2)
      real  :: ssacoice(jpb1:jpb2), forwice(jpb1:jpb2)
      real  :: extcoliq(jpb1:jpb2), gliq(jpb1:jpb2)
      real  :: ssacoliq(jpb1:jpb2), forwliq(jpb1:jpb2)

! Initialize

      hvrcld = '$Revision$'

!      ncbands = 29
      ib1 = jpb1
      ib2 = jpb2
      tauctot(:) = 0. 

      do lay = 1, nlayers
         do ib = ib1 , ib2
            taucldorig(lay,ib) = tauc(ib-15,lay)
            taucloud(lay,ib) = 0.0 
            ssacloud(lay,ib) = 1.0 
            asmcloud(lay,ib) = 0.0 
            tauctot(lay) = tauctot(lay) + tauc(ib-15,lay)
         enddo
      enddo

! Main layer loop
      do lay = 1, nlayers

         cwp = ciwp(lay) + clwp(lay)
         if (cldfrac(lay) .ge. cldmin .and. &
            (cwp .ge. cldmin .or. tauctot(lay) .ge. cldmin)) then

! (inflag=0): Cloud optical properties input directly
! Cloud optical properties already defined in tauc, ssac, asmc are unscaled;
! Apply delta-M scaling here
            if (inflag .eq. 0) then

               do ib = ib1 , ib2
                  taucldorig_a = tauc(ib-15,lay)
                  ffp = fsfc(ib-15,lay)
                  ffp1 = 1.0  - ffp
                  ffpssa = 1.0  - ffp * ssac(ib-15,lay)
                  ssacloud_a = ffp1 * ssac(ib-15,lay) / ffpssa
                  taucloud_a = ffpssa * taucldorig_a

                  taucldorig(lay,ib) = taucldorig_a
                  ssacloud(lay,ib) = ssacloud_a
                  taucloud(lay,ib) = taucloud_a
                  asmcloud(lay,ib) = (asmc(ib-15,lay) - ffp) / (ffp1)
               enddo

! (inflag=2): Separate treatement of ice clouds and water clouds.
            elseif (inflag .eq. 2) then       
               radice = rei(lay)

! Calculation of absorption coefficients due to ice clouds.
               if (ciwp(lay) .eq. 0.0 ) then
                  do ib = ib1 , ib2
                     extcoice(ib) = 0.0 
                     ssacoice(ib) = 0.0 
                     gice(ib)     = 0.0 
                     forwice(ib)  = 0.0 
                  enddo

! (iceflag = 1): 
! Note: This option uses Ebert and Curry approach for all particle sizes similar to
! CAM3 implementation, though this is somewhat ineffective for large ice particles
               elseif (iceflag .eq. 1) then
                  if (radice .lt. 13.0  .or. radice .gt. 130. ) stop &
                     'ICE RADIUS OUT OF BOUNDS'
                  do ib = ib1, ib2
                     if (wavenum2(ib) .gt. 1.43e04 ) then
                        icx = 1
                     elseif (wavenum2(ib) .gt. 7.7e03 ) then
                        icx = 2
                     elseif (wavenum2(ib) .gt. 5.3e03 ) then
                        icx = 3
                     elseif (wavenum2(ib) .gt. 4.0e03 ) then
                        icx = 4
                     elseif (wavenum2(ib) .ge. 2.5e03 ) then
                        icx = 5
                     endif
                     extcoice(ib) = abari(icx) + bbari(icx)/radice
                     ssacoice(ib) = 1.  - cbari(icx) - dbari(icx) * radice
                     gice(ib) = ebari(icx) + fbari(icx) * radice

! Check to ensure upper limit of gice is within physical limits for large particles
                     if (gice(ib) .ge. 1.0 ) gice(ib) = 1.0  - eps
                     forwice(ib) = gice(ib)*gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoice(ib) .lt. 0.0 ) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ib) .gt. 1.0 ) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ib) .lt. 0.0 ) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ib) .gt. 1.0 ) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ib) .lt. 0.0 ) stop 'ICE ASYM LESS THAN 0.0'
                  enddo

! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns

               elseif (iceflag .eq. 2) then
                  if (radice .lt. 5.0  .or. radice .gt. 131.0 ) stop 'ICE RADIUS OUT OF BOUNDS'
                  factor = (radice - 2. )/3. 
                  index = int(factor)
                  if (index .eq. 43) index = 42
                  fint = factor - float(index)
                  do ib = ib1, ib2
                     extcoice(ib) = extice2(index,ib) + fint * &
                                   (extice2(index+1,ib) -  extice2(index,ib))
                     ssacoice(ib) = ssaice2(index,ib) + fint * &
                                   (ssaice2(index+1,ib) -  ssaice2(index,ib))
                     gice(ib) = asyice2(index,ib) + fint * &
                                   (asyice2(index+1,ib) -  asyice2(index,ib))
                     forwice(ib) = gice(ib)*gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoice(ib) .lt. 0.0 ) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ib) .gt. 1.0 ) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ib) .lt. 0.0 ) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ib) .gt. 1.0 ) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ib) .lt. 0.0 ) stop 'ICE ASYM LESS THAN 0.0'
                  enddo

! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

               elseif (iceflag .eq. 3) then
                  if (radice .lt. 5.0  .or. radice .gt. 140.0 ) stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
                  factor = (radice - 2. )/3. 
                  index = int(factor)
                  if (index .eq. 46) index = 45
                  fint = factor - float(index)
                  do ib = ib1 , ib2
                     extcoice(ib) = extice3(index,ib) + fint * &
                                   (extice3(index+1,ib) - extice3(index,ib))
                     ssacoice(ib) = ssaice3(index,ib) + fint * &
                                   (ssaice3(index+1,ib) - ssaice3(index,ib))
                     gice(ib) = asyice3(index,ib) + fint * &
                               (asyice3(index+1,ib) - asyice3(index,ib))
                     fdelta(ib) = fdlice3(index,ib) + fint * &
                                 (fdlice3(index+1,ib) - fdlice3(index,ib))
                     if (fdelta(ib) .lt. 0.0 ) stop 'FDELTA LESS THAN 0.0'
                     if (fdelta(ib) .gt. 1.0 ) stop 'FDELTA GT THAN 1.0'                     
                     forwice(ib) = fdelta(ib) + 0.5  / ssacoice(ib)
! See Fu 1996 p. 2067 
                     if (forwice(ib) .gt. gice(ib)) forwice(ib) = gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoice(ib) .lt. 0.0 ) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ib) .gt. 1.0 ) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ib) .lt. 0.0 ) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ib) .gt. 1.0 ) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ib) .lt. 0.0 ) stop 'ICE ASYM LESS THAN 0.0'
                  enddo

               endif
                  
! Calculation of absorption coefficients due to water clouds.
                if (clwp(lay) .eq. 0.0 ) then
                   do ib = ib1 , ib2
                      extcoliq(ib) = 0.0 
                      ssacoliq(ib) = 0.0 
                      gliq(ib) = 0.0 
                      forwliq(ib) = 0.0 
                   enddo

                elseif (liqflag .eq. 1) then
                   radliq = rel(lay)
                   if (radliq .lt. 2.5  .or. radliq .gt. 60. ) stop &
                      'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
                   index = int(radliq - 1.5 )
                   if (index .eq. 0) index = 1
                   if (index .eq. 58) index = 57
                   fint = radliq - 1.5  - float(index)
                   do ib = ib1 , ib2
                      extcoliq(ib) = extliq1(index,ib) + fint * &
                                    (extliq1(index+1,ib) - extliq1(index,ib))
                      ssacoliq(ib) = ssaliq1(index,ib) + fint * &
                                    (ssaliq1(index+1,ib) - ssaliq1(index,ib))
                      if (fint .lt. 0.  .and. ssacoliq(ib) .gt. 1. ) &
                                     ssacoliq(ib) = ssaliq1(index,ib)
                      gliq(ib) = asyliq1(index,ib) + fint * &
                                (asyliq1(index+1,ib) - asyliq1(index,ib))
                      forwliq(ib) = gliq(ib)*gliq(ib)
! Check to ensure all calculated quantities are within physical limits.
                      if (extcoliq(ib) .lt. 0.0 ) stop 'LIQUID EXTINCTION LESS THAN 0.0'
                      if (ssacoliq(ib) .gt. 1.0 ) stop 'LIQUID SSA GRTR THAN 1.0'
                      if (ssacoliq(ib) .lt. 0.0 ) stop 'LIQUID SSA LESS THAN 0.0'
                      if (gliq(ib) .gt. 1.0 ) stop 'LIQUID ASYM GRTR THAN 1.0'
                      if (gliq(ib) .lt. 0.0 ) stop 'LIQUID ASYM LESS THAN 0.0'
                   enddo
                endif

                do ib = ib1 , ib2
                   tauliqorig = clwp(lay) * extcoliq(ib)
                   tauiceorig = ciwp(lay) * extcoice(ib)
                   taucldorig(lay,ib) = tauliqorig + tauiceorig

                   ssaliq = ssacoliq(ib) * (1.0  - forwliq(ib)) / &
                           (1.0  - forwliq(ib) * ssacoliq(ib))
                   tauliq = (1.0  - forwliq(ib) * ssacoliq(ib)) * tauliqorig
                   ssaice = ssacoice(ib) * (1.0  - forwice(ib)) / &
                           (1.0  - forwice(ib) * ssacoice(ib))
                   tauice = (1.0  - forwice(ib) * ssacoice(ib)) * tauiceorig

                   scatliq = ssaliq * tauliq
                   scatice = ssaice * tauice

                   taucloud(lay,ib) = tauliq + tauice

! Ensure non-zero taucmc and scatice
                   if (taucloud(lay,ib).eq.0.0 ) taucloud(lay,ib) = cldmin
                   if (scatice.eq.0.0 ) scatice = cldmin

                   ssacloud(lay,ib) = (scatliq + scatice) / taucloud(lay,ib)

                   if (iceflag .eq. 3) then
! In accordance with the 1996 Fu paper, equation A.3, 
! the moments for ice were calculated depending on whether using spheres
! or hexagonal ice crystals.
                      istr = 1
                      asmcloud(lay,ib) = (1.0 /(scatliq+scatice)) * &
                         (scatliq*(gliq(ib)**istr - forwliq(ib)) / &
                         (1.0  - forwliq(ib)) + scatice * ((gice(ib)-forwice(ib)) / &
                         (1.0  - forwice(ib)))**istr)
                   else 
! This code is the standard method for delta-m scaling. 
                      istr = 1
                      asmcloud(lay,ib) = (scatliq *  &
                         (gliq(ib)**istr - forwliq(ib)) / &
                         (1.0  - forwliq(ib)) + scatice * (gice(ib)**istr - forwice(ib)) / &
                         (1.0  - forwice(ib)))/(scatliq + scatice)
                   endif 

                enddo

            endif

         endif

! End layer loop
      enddo

      end subroutine cldprop_sw

      end module rrtmg_sw_cldprop


