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
   use rrlw_wvn, only: ngb
   use rrlw_cld, only: &
     absliq1, absice0, absice1, absice2, absice3, absice4, &
     ice1b

   implicit none

contains
	  
   ! ------------------------------------------------------------------------------
   subroutine cldprmc (ncol, nlay, &
                      cldfmc, ciwpmc, clwpmc, reice, reliq, &
                      iceflag, liqflag, taucmc, cloudy)
   ! ------------------------------------------------------------------------------
   ! Compute the cloud optical depths for each cloudy layer
   ! The g-point indices are really McICA subcolumns

      integer, intent(in) :: ncol  ! number of columns
      integer, intent(in) :: nlay  ! number of layers
      
      real,    intent(in) :: cldfmc (nlay,ngptlw,ncol)    ! cloud fraction [mcica]
      real,    intent(in) :: ciwpmc (nlay,ngptlw,ncol)    ! in-cloud ice water path [mcica]
      real,    intent(in) :: clwpmc (nlay,ngptlw,ncol)    ! in-cloud liq water path [mcica]

      real,    intent(in) :: reice (nlay,ncol)            ! ice crystal effective size [um]
      real,    intent(in) :: reliq (nlay,ncol)            ! liq droplet effective radius [um]

      integer, intent(in) :: iceflag, liqflag             ! cloud optical depth methods

      real,    intent(out) :: taucmc (nlay,ngptlw,ncol)   ! cloud optical depth [mcica]
      logical, intent(out) :: cloudy (nlay,ncol)          ! gridbox cloudy for any gpoint

      ! ------- Local -------

      integer :: icol      ! gridcolumn index
      integer :: ilay      ! layer index
      integer :: ig        ! g-point interval index
      integer :: ib        ! spectral band index

      real    :: radice    ! cloud ice effective size (microns)
      real    :: abscoice  ! ice absorption coefficients
      real    :: abscoliq  ! liquid absorption coefficients

      ! table lookup parameters
      integer :: index
      real :: factor, fint

      ! tiny threshold value for cloud water path, at or below
      ! which ice or liquid are treated as optically inactive.
      real, parameter :: cwp_negligible = 1.e-20

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
!                    Liquid effective radius limited to [2.5,60.0] um.

      ! flag gridbox if cloudy for any g-point
      cloudy = .false.
      do icol = 1,ncol
         do ilay = 1,nlay
            do ig = 1,ngptlw
               if (cldfmc(ilay,ig,icol) > 0.) then
                  cloudy(ilay,icol) = .true.
                  exit
               endif
            end do
         end do
      end do

      ! zero cloud optical thickness default
      taucmc = 0.

      ! Calculation of absorption coefficients due to ice clouds
      if (iceflag == 0) then

         do icol = 1,ncol
            do ilay = 1,nlay
               if (cloudy(ilay,icol)) then

                  ! 1 cloud band
                  abscoice = absice0(1) + absice0(2) / reice(ilay,icol)

                  do ig = 1,ngptlw
                     if (cldfmc(ilay,ig,icol) > 0. &
                        .and. ciwpmc(ilay,ig,icol) > cwp_negligible) then

                        taucmc(ilay,ig,icol) = ciwpmc(ilay,ig,icol) * abscoice

                     end if
                  end do  ! g-point

               end if  ! cloudy for any gpoint
            end do  ! layer
         end do  ! column

      elseif (iceflag == 1) then

         do icol = 1,ncol
            do ilay = 1,nlay
               if (cloudy(ilay,icol)) then

                  radice = reice(ilay,icol)

                  do ig = 1,ngptlw
                     if (cldfmc(ilay,ig,icol) > 0. &
                        .and. ciwpmc(ilay,ig,icol) > cwp_negligible) then

                        ! 5 cloud bands
                        ib = ice1b(ngb(ig))
                        abscoice = absice1(1,ib) + absice1(2,ib) / radice
                        taucmc(ilay,ig,icol) = ciwpmc(ilay,ig,icol) * abscoice

                     endif
                  end do  ! g-point

               end if  ! cloudy for any gpoint
            end do  ! layer
         end do  ! column

      elseif (iceflag == 2) then

         do icol = 1,ncol
            do ilay = 1,nlay
               if (cloudy(ilay,icol)) then

                  ! ice particle effective radius [5.0,131.0] um
                  factor = (reice(ilay,icol) - 2.)/3. 
                  index = int(factor)
                  if (index >= 43) then
                     if (index == 43) then 
                        index = 42
                     else
                        error stop 'cldprmc: iceflag 2: excessive high-radius extrapolation forbidden!'
                     end if
                  else if (index <= 0) then
                     if (index == 0) then 
                        index = 1
                     else
                        error stop 'cldprmc: iceflag 2: excessive low-radius extrapolation forbidden!'
                     end if
                  end if
                  fint = factor - float(index)

                  do ig = 1,ngptlw
                     if (cldfmc(ilay,ig,icol) > 0. &
                        .and. ciwpmc(ilay,ig,icol) > cwp_negligible) then

                        ! 16 cloud bands
                        ib = ngb(ig)
                        abscoice = &
                           absice2(index,ib) + fint * &
                           (absice2(index+1,ib) - (absice2(index,ib))) 
                        taucmc(ilay,ig,icol) = ciwpmc(ilay,ig,icol) * abscoice

                     endif
                  end do  ! g-point

               end if  ! cloudy for any gpoint
            end do  ! layer
         end do  ! column

      elseif (iceflag == 3) then

         do icol = 1,ncol
            do ilay = 1,nlay
               if (cloudy(ilay,icol)) then

                  ! ice particle generalized effective size [5.0,140.0] um
                  factor = (reice(ilay,icol) - 2.)/3. 
                  index = int(factor)
                  if (index >= 46) then
                     if (index == 46) then 
                        index = 45
                     else
                        error stop 'cldprmc: iceflag 3: excessive high-radius extrapolation forbidden!'
                     end if
                  else if (index <= 0) then
                     if (index == 0) then 
                        index = 1
                     else
                        error stop 'cldprmc: iceflag 3: excessive low-radius extrapolation forbidden!'
                     end if
                  end if
                  fint = factor - float(index)

                  do ig = 1,ngptlw
                     if (cldfmc(ilay,ig,icol) > 0. &
                        .and. ciwpmc(ilay,ig,icol) > cwp_negligible) then

                        ! 16 cloud bands
                        ib = ngb(ig)
                        abscoice = &
                           absice3(index,ib) + fint * &
                           (absice3(index+1,ib) - (absice3(index,ib)))
                        taucmc(ilay,ig,icol) = ciwpmc(ilay,ig,icol) * abscoice

                     endif
                  end do  ! g-point

               end if  ! cloudy for any gpoint
            end do  ! layer
         end do  ! column

      elseif (iceflag == 4) then

         do icol = 1,ncol
            do ilay = 1,nlay
               if (cloudy(ilay,icol)) then

                  ! ice particle effective diameter [1.0,200.0] um
                  factor = reice(ilay,icol)
                  index = int(factor)
                  if (index >= 200) then
                     if (index == 200) then 
                        index = 199
                     else
                        error stop 'cldprmc: iceflag 4: excessive high-radius extrapolation forbidden!'
                     end if
                  else if (index <= 0) then
                     if (index == 0) then 
                        index = 1
                     else
                        error stop 'cldprmc: iceflag 4: excessive low-radius extrapolation forbidden!'
                     end if
                  end if
                  fint = factor - float(index)

                  do ig = 1,ngptlw
                     if (cldfmc(ilay,ig,icol) > 0. &
                        .and. ciwpmc(ilay,ig,icol) > cwp_negligible) then

                        ! 16 cloud bands
                        ib = ngb(ig)
                        abscoice = &
                           absice4(index,ib) + fint * &
                           (absice4(index+1,ib) - (absice4(index,ib)))
                        taucmc(ilay,ig,icol) = ciwpmc(ilay,ig,icol) * abscoice

                     endif
                  end do  ! g-point

               end if  ! cloudy for any gpoint
            end do  ! layer
         end do  ! column
                  
      else
         error stop 'cldprmc: invalid iceflag'
      endif
        
      ! Calculation of absorption coefficients due to water clouds
      if (liqflag == 1) then

         do icol = 1,ncol
            do ilay = 1,nlay
               if (cloudy(ilay,icol)) then

                  ! liquid effective radius [2.5,60.0] um
                  factor = reliq(ilay,icol) - 1.5
                  index = int(factor)
                  if (index >= 58) then
                     if (index == 58) then 
                        index = 57
                     else
                        error stop 'cldprmc: liqflag 1: excessive high-radius extrapolation forbidden!'
                     end if
                  else if (index <= 0) then
                     if (index == 0) then 
                        index = 1
                     else
                        error stop 'cldprmc: liqflag 1: excessive low-radius extrapolation forbidden!'
                     end if
                  end if
                  fint = factor - float(index)

                  do ig = 1,ngptlw
                     if (cldfmc(ilay,ig,icol) > 0. &
                           .and. clwpmc(ilay,ig,icol) > cwp_negligible) then

                        ! 16 cloud bands
                        ib = ngb(ig)
                        abscoliq = &
                           absliq1(index,ib) + fint * &
                           (absliq1(index+1,ib) - (absliq1(index,ib)))
                        taucmc(ilay,ig,icol) = taucmc(ilay,ig,icol) + &
                           clwpmc(ilay,ig,icol) * abscoliq

                     endif
                  end do  ! g-point

               end if  ! cloudy for any gpoint
            end do  ! layer
         end do  ! column

      else
         error stop 'cldprmc: invalid liqflag'
      endif

      ! refine cloudy flag (for greater speed in rtrnmc) so that
      ! flags gridbox if OPTICALLY cloudy for any g-point. Issue
      ! is that cloudy() layers defined above may still not have
      ! any OPTICALLY significantly cloud (all taucmc may be zero),
      ! in which case they should be reclassed to .not.cloudy().
      do icol = 1,ncol
         layer: do ilay = 1,nlay
            if (cloudy(ilay,icol)) then 
               do ig = 1,ngptlw
                  ! if at least one g-point is optically cloudy
                  !    then layer remains cloudy
                  if (taucmc(ilay,ig,icol) > 0.) cycle layer
               end do
               ! no optically cloudy g-points so reclassed not cloudy
               cloudy(ilay,icol) = .false. 
            end if
         end do layer
      end do

   end subroutine cldprmc

end module rrtmg_lw_cldprmc
