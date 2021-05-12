! (dmb 2012) This is the GPU version of the cldprmc routine.  I have parallelized across 
! all 3 dimensions (columns, g-points, and layers) to make this routine run very fast on the GPU.  
! The greatest speedup was obtained by switching the indices for the cloud variables so that 
! the columns were the least significant (leftmost) dimension

module rrtmg_lw_cldprmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

   use parrrtm, only : ngptlw
   use rrlw_cld, only: &
     absliq1, absice0, absice1, absice2, absice3, absice4, &
     ice1b

   implicit none

contains
	  
   ! ------------------------------------------------------------------------------
   subroutine cldprmc(ncol, nlay, &
                      cldfmc, ciwpmc, clwpmc, reice, reliq, &
                      iceflag, liqflag, ngb, &
                      taucmc, icldlyr)
   ! ------------------------------------------------------------------------------
   ! Compute the cloud optical depths for each cloudy layer
   ! The g-point indices are really McICA subcolumns

      integer, intent(in) :: ncol  ! number of columns
      integer, intent(in) :: nlay  ! number of layers
      
      real,    intent(in) :: cldfmc(ncol,ngptlw,nlay)  ! cloud fraction [mcica]
      real,    intent(in) :: ciwpmc(ncol,ngptlw,nlay)  ! in-cloud ice water path [mcica]
      real,    intent(in) :: clwpmc(ncol,ngptlw,nlay)  ! in-cloud liquid water path [mcica]

      real,    intent(in) :: reice(ncol,nlay)          ! ice crystal effective size [um]
      real,    intent(in) :: reliq(ncol,nlay)          ! liq droplet effective radius [um]

      integer, intent(in) :: iceflag, liqflag            ! cloud optical depth methods
      integer, intent(in) :: ngb(ngptlw)

      real,    intent(out) :: taucmc(ncol,ngptlw,nlay)   ! cloud optical depth [mcica]
      integer, intent(out) :: icldlyr(ncol,nlay)         ! cloudy in any layer of COLUMN

      ! ------- Local -------

      integer :: icol                      ! gridcolumn index
      integer :: ilay                      ! layer index
      integer :: ib                        ! spectral band index
      integer :: ig                        ! g-point interval index

      real    :: cwp                       ! cloud water path
      real    :: radice                    ! cloud ice effective size (microns)
      real    :: radliq                    ! cloud liquid droplet radius (microns)
      real    :: abscoice                  ! ice absorption coefficients
      real    :: abscoliq                  ! liquid absorption coefficients

      ! table lookup parameters
      integer :: index
      real :: factor, fint

      real, parameter :: cldmin = 1.e-20   ! minimum value for cloud fraction
      real, parameter :: cwpmin = 1.e-20   ! minimum value for cloud water path

! ------- Definitions -------

!     Optical properties method distinguishes between liquid and ice clouds,
!     and requires further user input to specify the method to be used to 
!     compute the aborption due to each. For each cloudy layer, the cloud
!     fraction, and cloud liquid and ice water paths (g/m2) are input.
!
!       ICEFLAG = 0: The ice effective radius (microns) is input and the
!                    optical depths due to ice clouds are computed as in CCM3.
!                    Ice effective radius, r_ec, (Ebert and Curry, 1992),
!                      r_ec must be >= 10.0 microns
!
!       ICEFLAG = 1: The ice effective radius (microns) is input and the
!                    optical depths due to ice clouds are computed as in 
!                    Ebert and Curry, JGR, 97, 3831-3836 (1992).  The 
!                    spectral regions in this work have been matched with
!                    the spectral bands in RRTM to as great an extent 
!                    as possible:  
!                    E&C 1      IB = 5      RRTM bands 9-16
!                    E&C 2      IB = 4      RRTM bands 6-8
!                    E&C 3      IB = 3      RRTM bands 3-5
!                    E&C 4      IB = 2      RRTM band 2
!                    E&C 5      IB = 1      RRTM band 1
!                    Ice effective radius, r_ec, (Ebert and Curry, 1992),
!                      r_ec range is limited to 13.0 to 130.0 microns
!
!       ICEFLAG = 2: The ice effective radius (microns) is input and the
!                    optical properties due to ice clouds are computed from
!                    the optical properties stored in the RT code,
!                    STREAMER v3.0 (Reference: Key. J., Streamer 
!                    User's Guide, Cooperative Institute for
!                    Meteorological Satellite Studies, 2001, 96 pp.).
!                    Valid range of values for re are between 5.0 and
!                    131.0 micron.
!                    Ice effective radius, r_k, (Key, Streamer Ref. Manual,
!                      1996) r_k range is limited to 5.0 to 131.0 microns.
!
!       ICEFLAG = 3: The ice generalized effective size (dge) is input
!                    and the optical properties, are calculated as in
!                    Q. Fu, J. Climate, (1998). Q. Fu provided high resolution
!                    tables which were appropriately averaged for the
!                    bands in RRTM_LW.  Linear interpolation is used to
!                    get the coefficients from the stored tables.
!                    Valid range of values for dge are between 5.0 and
!                    140.0 micron.
!                    Generalized effective size, dge, (Fu, 1996),
!                      dge range is limited to 5.0 to 140.0 microns
!                      [dge = 1.0315 * r_ec]
!
!       LIQFLAG = 1: The water droplet effective radius (microns) is input 
!                    and the optical depths due to water clouds are computed 
!                    as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
!                    The values for absorption coefficients appropriate for
!                    the spectral bands in RRTM have been obtained for a 
!                    range of effective radii by an averaging procedure 
!                    based on the work of J. Pinto (private communication).
!                    Linear interpolation is used to get the absorption 
!                    coefficients for the input effective radius.

      icldlyr = 0
      taucmc = 0.
      do icol = 1,ncol
         do ilay = 1,nlay
            do ig = 1,ngptlw

               ! (dmb 2012) all of the cloud variables have been modified
               !   so that the column dimensions is least significant.

               if (cldfmc(icol,ig,ilay) == 1. ) then
                  ! flag if layer is cloudy in any layer of the column
                  ! i.e., for any subcolumn (g-point) at that layer
                  icldlyr(icol,ilay) = 1
               endif
               cwp = ciwpmc(icol,ig,ilay) + clwpmc(icol,ig,ilay)
   
               if (cldfmc(icol,ig,ilay) > cldmin .and. cwp > cwpmin) then

                  radice = reice(icol,ilay)

                  ! Calculation of absorption coefficients due to ice clouds.
                  if (ciwpmc(icol,ig,ilay) == 0.) then

                     abscoice = 0. 
	 
                  elseif (iceflag == 0) then

                     ! 1 cloud band
                     abscoice = absice0(1) + absice0(2)/radice

                  elseif (iceflag == 1) then

                     ! 5 cloud bands
                     ib = ice1b(ngb(ig))
                     abscoice = absice1(1,ib) + absice1(2,ib)/radice

                  ! For iceflag=2 option, ice particle effective radius
                  !   is limited to 5.0 to 131.0 microns
                  elseif (iceflag == 2) then

                     ! 16 cloud bands
                     factor = (radice - 2.)/3. 
                     index = int(factor)
                     if (index == 43) index = 42
                     fint = factor - float(index)
                     ib = ngb(ig)
                     abscoice = &
                        absice2(index,ib) + fint * &
                        (absice2(index+1,ib) - (absice2(index,ib))) 
            
                  ! For iceflag=3 option, ice particle generalized effective
                  !   size is limited to 5.0 to 140.0 microns
                  elseif (iceflag == 3) then

                     ! 16 cloud bands
                     factor = (radice - 2.)/3. 
                     index = int(factor)
                     if (index == 46) index = 45
                     fint = factor - float(index)
                     ib = ngb(ig)
                     abscoice = &
                        absice3(index,ib) + fint * &
                        (absice3(index+1,ib) - (absice3(index,ib)))

                  ! For iceflag=4 option, ice particle effective diameter
                  !   is limited to 1.0 to 200.0 microns
                  elseif (iceflag == 4) then

                     ! 16 cloud bands
                     factor = radice
                     index = int(factor)
                     fint = factor - float(index)
                     ib = ngb(ig)
                     abscoice = &
                        absice4(index,ib) + fint * &
                        (absice4(index+1,ib) - (absice4(index,ib)))

                  else
                     error stop 'cldprmc: invalid iceflag'
                  endif
                  
                  ! Calculation of absorption coefficients due to water clouds.
                  if (clwpmc(icol,ig,ilay) == 0.) then

                     abscoliq = 0. 
	 
                  elseif (liqflag == 1) then

                     ! 16 cloud bands
                     radliq = reliq(icol,ilay)
                     index = int(radliq - 1.5)
                     if (index == 0) index = 1
                     if (index == 58) index = 57
                     fint = radliq - 1.5  - float(index)
                     ib = ngb(ig)
                     abscoliq = &
                        absliq1(index,ib) + fint * &
                        (absliq1(index+1,ib) - (absliq1(index,ib)))

                  else
                     error stop 'cldprmc: invalid liqflag'
                  endif

                  taucmc(icol,ig,ilay) = ciwpmc(icol,ig,ilay) * abscoice + &
                                         clwpmc(icol,ig,ilay) * abscoliq

               endif
        
            end do
         end do
      end do

   end subroutine cldprmc

end module rrtmg_lw_cldprmc
