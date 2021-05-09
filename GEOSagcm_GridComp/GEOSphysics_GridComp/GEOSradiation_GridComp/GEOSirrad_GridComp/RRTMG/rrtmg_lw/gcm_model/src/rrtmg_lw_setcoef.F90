!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

module rrtmg_lw_setcoef

   use parrrtm,  only : nbndlw
   use rrlw_wvn, only : totplnk, totplk16, totplnkderiv, totplk16deriv
   use rrlw_con, only : grav, avogad
   use rrlw_ref

   implicit none

   ! lower / upper atmosphere separator (see code)
   integer, allocatable :: laytrop (:)

   ! molecular concentrations (scaled molecules per cm^2)
   real, allocatable, dimension (:,:) :: &
      colh2o, colco2, colo3, coln2o, colch4, colo2, colco, colbrd, &
      colcfc11, colcfc12, colcfc22, colccl4

   real, allocatable :: coldry (:,:)  ! dry air column density [molec/cm^2]
   real, allocatable :: wbroad (:,:)  ! broadening gas column density [molec/cm^2]
   real, allocatable :: pwvcm  (:)    ! precipitable water vapor [cm]

   ! Planck functions and interpolation factors for continuum,
   ! minor absorber, and (ln p, T, & binary species) CKD tables
   real, allocatable, dimension (:,:,:) :: planklev, planklay
   real, allocatable, dimension (:,:) :: &
      plankbnd, dplankbnd_dt, &
      forfac, forfrac, selffac, selffrac, &
      scaleminor, scaleminorn2, minorfrac
   integer, allocatable, dimension (:,:) :: &
      jp, jt, jt1, indfor, indself, indminor
   real, allocatable, dimension (:,:) :: &
      rat_h2oco2, rat_h2oco2_1, rat_h2oo3,  rat_h2oo3_1,  &
      rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
      rat_n2oco2, rat_n2oco2_1, rat_o3co2,  rat_o3co2_1
   real, allocatable, dimension(:,:) :: &
      fac00, fac01, fac10, fac11

contains

   !----------------------------------------------------------------------------
   subroutine setcoef ( &
      ncol, nlay, istart, idrv, &
      pavel, tavel, pz, tz, tbound, semiss, &
      h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, covmr, &
      cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr)

   !----------------------------------------------------------------------------
   !
   ! For a given atmosphere, calculate the indices and fractions related to
   ! pressure and temperature interpolations. Also calculate the values of
   ! the integrated Planck functions for each band at the level and layer
   ! temperatures.

      integer, intent(in) :: ncol    ! number of gridcolumns
      integer, intent(in) :: nlay    ! number of layers
      integer, intent(in) :: istart  ! see code
      integer, intent(in) :: idrv    ! Planck derivative option flag

      real, intent(in) :: pavel  (ncol,nlay)    ! layer pressures [hPa]
      real, intent(in) :: tavel  (ncol,nlay)    ! layer temperatures [K]
      real, intent(in) :: pz     (ncol,0:nlay)  ! level (interface) pressure [hPa]
      real, intent(in) :: tz     (ncol,0:nlay)  ! level (interface) temperatures [K]
      real, intent(in) :: tbound (ncol)         ! surface temperature [K]
      real, intent(in) :: semiss (ncol,nbndlw)  ! surface emissivity

      ! volume mixing ratios (moles per mole of dry air)
      real, intent(in), dimension (ncol,nlay) :: &
         h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, covmr, &
         cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr

      real, parameter :: amd = 28.9660  ! Effective molecular weight of dry air (g/mol)
      real, parameter :: amw = 18.0160  ! Molecular weight of water vapor (g/mol)
      real, parameter :: stpfac = 296. / 1013.

      integer :: icol, lay, iband
      integer :: indbound, indlev0, indlay, indlev, jp1
      real    :: tbndfrac, t0frac, tlayfrac, tlevfrac
      real    :: dbdtlev, dbdtlay
      real    :: plog, fp, ft, ft1, water, scalefac, factor, compfp
      real    :: amm, amttl, wvttl, wvsh, summol, btemp, wv, lcoldry

      ! allocate or reallocate as necessary
      if (allocated (laytrop)) then
         if (size(laytrop) /= ncol) then
           deallocate (coldry      )
           deallocate (wbroad      )
           deallocate (pwvcm       )
           deallocate (laytrop     )
           deallocate (colh2o      )
           deallocate (colco2      )
           deallocate (colo3       )
           deallocate (coln2o      )
           deallocate (colch4      )
           deallocate (colo2       )
           deallocate (colco       )
           deallocate (colbrd      )
           deallocate (colcfc11    )
           deallocate (colcfc12    )
           deallocate (colcfc22    )
           deallocate (colccl4     )
           deallocate (plankbnd    )
           deallocate (dplankbnd_dt)
           deallocate (planklev    )
           deallocate (planklay    )
           deallocate (jp          )
           deallocate (jt          )
           deallocate (jt1         )
           deallocate (forfac      )
           deallocate (indfor      )
           deallocate (forfrac     )
           deallocate (selffac     )
           deallocate (indself     )
           deallocate (selffrac    )
           deallocate (scaleminor  )
           deallocate (scaleminorn2)
           deallocate (indminor    )
           deallocate (minorfrac   )
           deallocate (rat_h2oco2  )
           deallocate (rat_h2oco2_1)
           deallocate (rat_h2oo3   )
           deallocate (rat_h2oo3_1 )
           deallocate (rat_h2on2o  )
           deallocate (rat_h2on2o_1)
           deallocate (rat_h2och4  )
           deallocate (rat_h2och4_1)
           deallocate (rat_n2oco2  )
           deallocate (rat_n2oco2_1)
           deallocate (rat_o3co2   )
           deallocate (rat_o3co2_1 )
           deallocate (fac00       )
           deallocate (fac01       )
           deallocate (fac10       )
           deallocate (fac11       )
         end if
      end if
      if (.not. allocated (laytrop)) then
         allocate (coldry       (ncol,nlay))
         allocate (wbroad       (ncol,nlay))
         allocate (pwvcm        (ncol))
         allocate (laytrop      (ncol))
         allocate (colh2o       (ncol,nlay))
         allocate (colco2       (ncol,nlay))
         allocate (colo3        (ncol,nlay))
         allocate (coln2o       (ncol,nlay))
         allocate (colch4       (ncol,nlay))
         allocate (colo2        (ncol,nlay))
         allocate (colco        (ncol,nlay))
         allocate (colbrd       (ncol,nlay))
         allocate (colcfc11     (ncol,nlay))
         allocate (colcfc12     (ncol,nlay))
         allocate (colcfc22     (ncol,nlay))
         allocate (colccl4      (ncol,nlay))
         allocate (plankbnd     (ncol,nbndlw))
         allocate (dplankbnd_dt (ncol,nbndlw))
         allocate (planklev     (ncol,0:nlay,nbndlw))
         allocate (planklay     (ncol,  nlay,nbndlw))
         allocate (jp           (ncol,nlay))
         allocate (jt           (ncol,nlay))
         allocate (jt1          (ncol,nlay))
         allocate (forfac       (ncol,nlay))
         allocate (indfor       (ncol,nlay))
         allocate (forfrac      (ncol,nlay))
         allocate (selffac      (ncol,nlay))
         allocate (indself      (ncol,nlay))
         allocate (selffrac     (ncol,nlay))
         allocate (scaleminor   (ncol,nlay))
         allocate (scaleminorn2 (ncol,nlay))
         allocate (indminor     (ncol,nlay))
         allocate (minorfrac    (ncol,nlay))
         allocate (rat_h2oco2   (ncol,nlay))
         allocate (rat_h2oco2_1 (ncol,nlay))
         allocate (rat_h2oo3    (ncol,nlay))
         allocate (rat_h2oo3_1  (ncol,nlay))
         allocate (rat_h2on2o   (ncol,nlay))
         allocate (rat_h2on2o_1 (ncol,nlay))
         allocate (rat_h2och4   (ncol,nlay))
         allocate (rat_h2och4_1 (ncol,nlay))
         allocate (rat_n2oco2   (ncol,nlay))
         allocate (rat_n2oco2_1 (ncol,nlay))
         allocate (rat_o3co2    (ncol,nlay))
         allocate (rat_o3co2_1  (ncol,nlay))
         allocate (fac00        (ncol,nlay))
         allocate (fac01        (ncol,nlay))
         allocate (fac10        (ncol,nlay))
         allocate (fac11        (ncol,nlay))
      end if

      ! (re-)initialize as necessary
      ! (if not done below for all entries)
      
      ! set coldry, wbroad, and pwvcm

      do icol = 1,ncol

         do lay = 1, nlay

            ! molecular weight of moist air in g/mol:
            !   amm = [(\sum_i n_i a_i) + n_w a_w] / [(\sum_i n_i) + n_w],
            ! where i is a dry air species, n_i and a_i are its number of
            ! moles (in a volume of air) nand molecular weight, and similarly
            ! n_w and a_w (== amw) for water molecules. In the same way,
            !   amd = [\sum_i n_i a_i] / n_d, where n_d == [\sum_i n_i],
            ! so,
            !   amm = [n_d amd + n_w amw] / [n_d + n_w]
            !       = [amd + h2ovmr amw] / [1 + h2ovmr]
            ! since h2ovmr = n_w / n_d.
            ! Finally, since h2ovmr ~ 1%, on average, at sea-level, less above,
            ! a good *approximation* is
            !   amm \approx [1 - h2ovmr] * amd + h2ovmr * amw.

            amm = (1. - h2ovmr(icol,lay)) * amd + h2ovmr(icol,lay) * amw

            ! coldry = layer molecules of dry air per cm^2:
            ! -dp/g = rho*dz, where -dp = 100*(pz(k-1)-pz(k)) since pz in [hPa].
            ! So 1e3*(pz(k-1)-pz(k))/(1e2*grav) = (g moist air)/cm^2. This /amm
            ! is (moles moist air)/cm^2, and *avogad is (molecules moist air)/cm^2.
            ! Finally number of dry air (non-water) molecules per cm^2 should be
            ! this times n_d / (n_d + n_w) = 1 / (1 + h2ovmr).

            coldry(icol,lay) = (pz(icol,lay-1)-pz(icol,lay)) * 1.e3 * avogad / &
                  (1.e2 * grav * amm * (1. + h2ovmr(icol,lay)))

         end do

         ! calculate precipitable water vapor, and
         ! calculate n2 column density (for broadening)

         amttl = 0.
         wvttl = 0.
         do lay = 1,nlay

            ! pmn: summol appears to be used as a proxy for the sum of
            ! non-n2 dry air mole fractions. In fact it also contains
            ! the non-radiatively active Argon, about 1% by volume.
            ! pmn: Does this need to be corrected?
            summol = co2vmr(icol,lay) + o3vmr(icol,lay) &
               + n2ovmr(icol,lay) + ch4vmr(icol,lay) + o2vmr(icol,lay)
            ! wbroad is ~ number of n2 molecules per cm^2 in layer
            wbroad(icol,lay) = coldry(icol,lay) * (1. - summol)

            ! btemp is number of h2o molecules per cm^2 in layer
            ! so wvttl accumulates column h2o molecules per cm^2
            ! and amttl the total column molecules per cm^2
            btemp = h2ovmr(icol,lay) * coldry(icol,lay)
            amttl = amttl + coldry(icol,lay) + btemp
            wvttl = wvttl + btemp

         enddo

         ! column specic humidity [grams water vapor / grams moist air]
         ! pmn: almost. Really needs amm not amd in denominator.
         wvsh = (amw * wvttl) / (amd * amttl)

         ! precipitable water vapor [grams water / cm^2], or,
         ! since 1 cm^3 of liquid water ~ 1 gram, it is also the
         ! precipitable water vapor in [cm].
         pwvcm(icol) = wvsh * (1.e3 * pz(icol,0)) / (1.e2 * grav)

      enddo

      ! main loop over column
      do icol = 1,ncol

         ! tbound indexing within Planck table
         indbound = tbound(icol) - 159. 
         if (indbound < 1) then
            indbound = 1
         elseif (indbound > 180) then
            indbound = 180
         endif
         tbndfrac = tbound(icol) - 159. - float(indbound)

         ! tz(lev==0) indexing within Planck table
         indlev0 = tz(icol,0) - 159. 
         if (indlev0 < 1) then
            indlev0 = 1
         elseif (indlev0 > 180) then
            indlev0 = 180
         endif
         t0frac = tz(icol,0) - 159. - float(indlev0)

         ! laytrop is the highest layer with pavel > ~95.6 hPa
         ! (see below for details)
         laytrop(icol) = 0

         ! layer loop (bottom up)
         do lay = 1,nlay

            lcoldry = coldry(icol,lay)       ! dry air molec / cm^2
            wv = h2ovmr(icol,lay) * lcoldry  ! water vapor molec / cm^2

            ! layer temperature indexing within Planck table
            indlay = tavel(icol,lay) - 159. 
            if (indlay < 1) then
               indlay = 1
            elseif (indlay > 180) then
               indlay = 180
            endif
            tlayfrac = tavel(icol,lay) - 159. - float(indlay)

            ! level (interface) temperature indexing within Planck table
            indlev = tz(icol,lay) - 159. 
            if (indlev < 1) then
               indlev = 1
            elseif (indlev > 180) then
               indlev = 180
            endif
            tlevfrac = tz(icol,lay) - 159. - float(indlev)

            ! Calculate the integrated Planck functions for each band at the
            ! surface, level, and layer temperatures.

            ! spectral band loop 
            do iband = 1,15

               ! interpolate Planck functions to boundary & level zero temperatures
               if (lay.eq.1) then
                  dbdtlev = totplnk(indbound+1,iband) - totplnk(indbound,iband)
                  plankbnd(icol,iband) = semiss(icol,iband) * &
                     (totplnk(indbound,iband) + tbndfrac * dbdtlev)
                  dbdtlev = totplnk(indlev0+1,iband)-totplnk(indlev0,iband)
                  planklev(icol,0,iband) = totplnk(indlev0,iband) + t0frac * dbdtlev
                  if (idrv .eq. 1) then 
                     dbdtlev = totplnkderiv(indbound+1,iband) - totplnkderiv(indbound,iband)
                     dplankbnd_dt(icol,iband) = semiss(icol,iband) * &
                        (totplnkderiv(indbound,iband) + tbndfrac * dbdtlev)
                  endif
               endif

               ! interpolate Planck functions to layer & level temperatures
               dbdtlev = totplnk(indlev+1,iband) - totplnk(indlev,iband)
               planklev(icol,lay,iband) = totplnk(indlev,iband) + tlevfrac * dbdtlev
               dbdtlay = totplnk(indlay+1,iband) - totplnk(indlay,iband)
               planklay(icol,lay,iband) = totplnk(indlay,iband) + tlayfrac * dbdtlay

            enddo

            ! For band 16, if radiative transfer will be performed on just
            ! this band, use integrated Planck values up to 3250 cm-1.  
            ! If radiative transfer will be performed across all 16 bands,
            ! then include in the integrated Planck values for this band
            ! contributions from 2600 cm-1 to infinity.
            iband = 16
            if (istart == 16) then
               if (lay == 1) then
                  dbdtlev = totplk16(indbound+1) - totplk16(indbound)
                  plankbnd(icol,iband) = semiss(icol,iband) * &
                     (totplk16(indbound) + tbndfrac * dbdtlev)
! pmn: 2021-05-05
! istart==16 has not been used so far. Just as well, because the following
! line was in error by using totplnk instead of totplk16
                  dbdtlev = totplk16(indlev0+1) - totplk16(indlev0)
                  planklev(icol,0,iband) = totplk16(indlev0) + t0frac * dbdtlev
                  if (idrv == 1) then
                     dbdtlev = totplk16deriv(indbound+1) - totplk16deriv(indbound)
                     dplankbnd_dt(icol,iband) = semiss(icol,iband) * &
                        (totplk16deriv(indbound) + tbndfrac * dbdtlev)
                  endif
               endif
               dbdtlev = totplk16(indlev+1) - totplk16(indlev)
               planklev(icol,lay,iband) = totplk16(indlev) + tlevfrac * dbdtlev
               dbdtlay = totplk16(indlay+1) - totplk16(indlay)
               planklay(icol,lay,iband) = totplk16(indlay) + tlayfrac * dbdtlay
            else
               if (lay == 1) then
                  dbdtlev = totplnk(indbound+1,iband) - totplnk(indbound,iband)
                  plankbnd(icol,iband) = semiss(icol,iband) * &
                     (totplnk(indbound,iband) + tbndfrac * dbdtlev)
                  dbdtlev = totplnk(indlev0+1,iband) - totplnk(indlev0,iband)
                  planklev(icol,0,iband) = totplnk(indlev0,iband) + t0frac * dbdtlev
                  if (idrv == 1) then 
                     dbdtlev = totplnkderiv(indbound+1,iband) - totplnkderiv(indbound,iband)
                     dplankbnd_dt(icol,iband) = semiss(icol,iband) * &
                        (totplnkderiv(indbound,iband) + tbndfrac * dbdtlev)
                  endif
               endif
               dbdtlev = totplnk(indlev+1,iband) - totplnk(indlev,iband)
               planklev(icol,lay,iband) = totplnk(indlev,iband) + tlevfrac * dbdtlev
               dbdtlay = totplnk(indlay+1,iband) - totplnk(indlay,iband)
               planklay(icol,lay,iband) = totplnk(indlay,iband) + tlayfrac * dbdtlay
            endif

            ! Find the two reference pressures on either side of the
            ! layer pressure. Store them in JP and JP1. Store in FP the
            ! fraction of the difference (in ln(pressure)) between these
            ! two values that the layer pressure lies.

            plog = alog(pavel(icol,lay))
            jp(icol,lay) = int(36. - 5*(plog+0.04))
            if (jp(icol,lay) < 1) then
               jp(icol,lay) = 1
            elseif (jp(icol,lay) > 58) then
               jp(icol,lay) = 58
            endif
            jp1 = jp(icol,lay) + 1
            fp = 5. *(preflog(jp(icol,lay)) - plog)

            ! Determine, for each reference pressure (JP and JP1), which
            ! reference temperature (these are different for each  
            ! reference pressure) is nearest the layer temperature but does
            ! not exceed it. Store these indices in JT and JT1, resp.
            ! Store in FT (resp. FT1) the fraction of the way between JT
            ! (JT1) and the next highest reference temperature that the 
            ! layer temperature falls.

            jt(icol,lay) = int(3. + (tavel(icol,lay)-tref(jp(icol,lay)))/15.)
            if (jt(icol,lay) < 1) then
               jt(icol,lay) = 1
            elseif (jt(icol,lay) > 4) then
               jt(icol,lay) = 4
            endif
            ft = ((tavel(icol,lay)-tref(jp(icol,lay)))/15.) - float(jt(icol,lay)-3)

            jt1(icol,lay) = int(3. + (tavel(icol,lay)-tref(jp1))/15.)
            if (jt1(icol,lay) < 1) then
               jt1(icol,lay) = 1
            elseif (jt1(icol,lay) > 4) then
               jt1(icol,lay) = 4
            endif
            ft1 = ((tavel(icol,lay)-tref(jp1))/15.) - float(jt1(icol,lay)-3)

            ! intermediates needed for selffac and forfac
            water = wv / lcoldry                                   ! same as h2ovmr
            scalefac = pavel(icol,lay) * stpfac / tavel(icol,lay)  ! normalized p/T scalefac

            ! If the pressure is less than ~100mb, perform a different
            ! set of species interpolations.

            ! lower atmosphere
            if (plog > 4.56) then

               ! So laytrop counts the layers with pavel > exp(4.56) ~ 95.6 hPa.
               ! But pavel is strictly decreasing with height in model atmosphere.
               ! So pavel(lay <= laytrop) > ~95.6 hPa and pavel(lay > laytrop) <= ~95.6 hPa
               laytrop(icol) = laytrop(icol) + 1

               ! factors needed for water vapor foreign-continuum
               forfac(icol,lay) = scalefac/(1.+water)
               factor = (332.-tavel(icol,lay))/36. 
               indfor(icol,lay) = min(2,max(1,int(factor)))
               forfrac(icol,lay) = factor - float(indfor(icol,lay))

               ! factors needed to separately include the water vapor
               ! self-continuum in the calc'n of absorption coefficient
               selffac(icol,lay) = water * forfac(icol,lay)
               factor = (tavel(icol,lay)-188.)/7.2 
               indself(icol,lay) = min(9,max(1,int(factor)-7))
               selffrac(icol,lay) = factor - float(indself(icol,lay)+7)

               ! factors needed to separately include the minor gases
               ! in the calculation of absorption coefficient
               scaleminor(icol,lay) = pavel(icol,lay)/tavel(icol,lay)
               scaleminorn2(icol,lay) = (pavel(icol,lay)/tavel(icol,lay)) &
                  * (wbroad(icol,lay)/(lcoldry+wv))
               factor = (tavel(icol,lay)-180.8)/7.2 
               indminor(icol,lay) = min(18,max(1,int(factor)))
               minorfrac(icol,lay) = factor - float(indminor(icol,lay))

               ! Setup reference ratio to be used in calculation of binary
               ! species parameter in lower atmosphere.
               rat_h2oco2  (icol,lay) = chi_mls(1,jp(icol,lay)  ) / chi_mls(2,jp(icol,lay)  )
               rat_h2oco2_1(icol,lay) = chi_mls(1,jp(icol,lay)+1) / chi_mls(2,jp(icol,lay)+1)

               rat_h2oo3   (icol,lay) = chi_mls(1,jp(icol,lay)  ) / chi_mls(3,jp(icol,lay)  )
               rat_h2oo3_1 (icol,lay) = chi_mls(1,jp(icol,lay)+1) / chi_mls(3,jp(icol,lay)+1)

               rat_h2on2o  (icol,lay) = chi_mls(1,jp(icol,lay)  ) / chi_mls(4,jp(icol,lay)  )
               rat_h2on2o_1(icol,lay) = chi_mls(1,jp(icol,lay)+1) / chi_mls(4,jp(icol,lay)+1)

               rat_h2och4  (icol,lay) = chi_mls(1,jp(icol,lay)  ) / chi_mls(6,jp(icol,lay)  )
               rat_h2och4_1(icol,lay) = chi_mls(1,jp(icol,lay)+1) / chi_mls(6,jp(icol,lay)+1)

               rat_n2oco2  (icol,lay) = chi_mls(4,jp(icol,lay)  ) / chi_mls(2,jp(icol,lay)  )
               rat_n2oco2_1(icol,lay) = chi_mls(4,jp(icol,lay)+1) / chi_mls(2,jp(icol,lay)+1)
!pmn: others must be initialized to zero? or to NaN if not needed

               ! Calculate needed column amounts.
               colh2o  (icol,lay) = 1.e-20 * h2ovmr  (icol,lay) * lcoldry
               colco2  (icol,lay) = 1.e-20 * co2vmr  (icol,lay) * lcoldry
               colo3   (icol,lay) = 1.e-20 * o3vmr   (icol,lay) * lcoldry
               coln2o  (icol,lay) = 1.e-20 * n2ovmr  (icol,lay) * lcoldry
               colch4  (icol,lay) = 1.e-20 * ch4vmr  (icol,lay) * lcoldry
               colo2   (icol,lay) = 1.e-20 * o2vmr   (icol,lay) * lcoldry
               colco   (icol,lay) = 1.e-20 * covmr   (icol,lay) * lcoldry
               colcfc11(icol,lay) = 1.e-20 * cfc11vmr(icol,lay) * lcoldry
               colcfc12(icol,lay) = 1.e-20 * cfc12vmr(icol,lay) * lcoldry
               colcfc22(icol,lay) = 1.e-20 * cfc22vmr(icol,lay) * lcoldry
               colccl4 (icol,lay) = 1.e-20 * ccl4vmr (icol,lay) * lcoldry
               colbrd  (icol,lay) = 1.e-20 * wbroad  (icol,lay)
               if (colco2(icol,lay) == 0.) colco2(icol,lay) = 1.e-32 * lcoldry
               if (colo3 (icol,lay) == 0.) colo3 (icol,lay) = 1.e-32 * lcoldry
               if (coln2o(icol,lay) == 0.) coln2o(icol,lay) = 1.e-32 * lcoldry
               if (colch4(icol,lay) == 0.) colch4(icol,lay) = 1.e-32 * lcoldry
               if (colco (icol,lay) == 0.) colco (icol,lay) = 1.e-32 * lcoldry

            else  ! upper atmosphere (above laytrop)

               ! factors needed for water vapor foreign-continuum
               forfac(icol,lay) = scalefac/(1.+water)
               factor = (tavel(icol,lay)-188.)/36. 
               indfor(icol,lay) = 3
               forfrac(icol,lay) = factor - 1.

               ! factors needed to separately include the water vapor
               ! self-continuum in the calc'n of absorption coefficient
               selffac(icol,lay) = water * forfac(icol,lay)
!pmn? others selfrac, indself
!pmn? perhaps self not needed in upper atmosphere, set to NaN?

               ! factors needed to separately include the minor gases
               ! in the calculation of absorption coefficient
               scaleminor(icol,lay) = pavel(icol,lay)/tavel(icol,lay)         
               scaleminorn2(icol,lay) = (pavel(icol,lay)/tavel(icol,lay)) &
                  * (wbroad(icol,lay)/(lcoldry+wv))
               factor = (tavel(icol,lay)-180.8)/7.2 
               indminor(icol,lay) = min(18,max(1,int(factor)))
               minorfrac(icol,lay) = factor - float(indminor(icol,lay))

               ! Setup reference ratio to be used in calculation of binary
               ! species parameter in upper atmosphere.
               rat_h2oco2  (icol,lay) = chi_mls(1,jp(icol,lay)  ) / chi_mls(2,jp(icol,lay)  )
               rat_h2oco2_1(icol,lay) = chi_mls(1,jp(icol,lay)+1) / chi_mls(2,jp(icol,lay)+1)

               rat_o3co2   (icol,lay) = chi_mls(3,jp(icol,lay)  ) / chi_mls(2,jp(icol,lay)  )
               rat_o3co2_1 (icol,lay) = chi_mls(3,jp(icol,lay)+1) / chi_mls(2,jp(icol,lay)+1)

! pmn? others must be initialized to zero? or NaN not needed in upper atmos?

               ! Calculate needed column amounts.
               colh2o  (icol,lay) = 1.e-20 * h2ovmr  (icol,lay) * lcoldry
               colco2  (icol,lay) = 1.e-20 * co2vmr  (icol,lay) * lcoldry
               colo3   (icol,lay) = 1.e-20 * o3vmr   (icol,lay) * lcoldry
               coln2o  (icol,lay) = 1.e-20 * n2ovmr  (icol,lay) * lcoldry
               colch4  (icol,lay) = 1.e-20 * ch4vmr  (icol,lay) * lcoldry
               colo2   (icol,lay) = 1.e-20 * o2vmr   (icol,lay) * lcoldry
               colco   (icol,lay) = 1.e-20 * covmr   (icol,lay) * lcoldry
               colcfc11(icol,lay) = 1.e-20 * cfc11vmr(icol,lay) * lcoldry
               colcfc12(icol,lay) = 1.e-20 * cfc12vmr(icol,lay) * lcoldry
               colcfc22(icol,lay) = 1.e-20 * cfc22vmr(icol,lay) * lcoldry
               colccl4 (icol,lay) = 1.e-20 * ccl4vmr (icol,lay) * lcoldry
               colbrd  (icol,lay) = 1.e-20 * wbroad  (icol,lay)
               if (colco2(icol,lay) == 0.) colco2(icol,lay) = 1.e-32 * lcoldry
               if (colo3 (icol,lay) == 0.) colo3 (icol,lay) = 1.e-32 * lcoldry
               if (coln2o(icol,lay) == 0.) coln2o(icol,lay) = 1.e-32 * lcoldry
               if (colch4(icol,lay) == 0.) colch4(icol,lay) = 1.e-32 * lcoldry
               if (colco (icol,lay) == 0.) colco (icol,lay) = 1.e-32 * lcoldry

            end if

            ! We have now isolated the layer ln pressure and temperature,
            ! between two reference pressures and two reference temperatures 
            ! (for each reference pressure).  We multiply the pressure 
            ! fraction FP with the appropriate temperature fractions to get 
            ! the factors that will be needed for the interpolation that yields
            ! the optical depths (performed in routines TAUGBn for band n).
            compfp = 1. - fp
            fac10(icol,lay) = compfp * ft
            fac00(icol,lay) = compfp * (1. - ft)
            fac11(icol,lay) = fp * ft1
            fac01(icol,lay) = fp * (1. - ft1)

            ! rescale selffac and forfac for use in taumol
            selffac(icol,lay) = colh2o(icol,lay) * selffac(icol,lay)
            forfac (icol,lay) = colh2o(icol,lay) * forfac (icol,lay)
	
         enddo ! layer
      end do ! column

   end subroutine setcoef


   ! Deallocate module level allocatables when no longer needed.
   ! Do this after only when rrtmg_lw() is done. It is not necessary
   ! to do this between partitions because setcoef() handles reallocation
   ! if the size of a partition changes.
   subroutine setcoef_free
      if (allocated (laytrop)) then
         deallocate (coldry      )
         deallocate (wbroad      )
         deallocate (pwvcm       )
         deallocate (laytrop     )
         deallocate (colh2o      )
         deallocate (colco2      )
         deallocate (colo3       )
         deallocate (coln2o      )
         deallocate (colch4      )
         deallocate (colo2       )
         deallocate (colco       )
         deallocate (colbrd      )
         deallocate (colcfc11    )
         deallocate (colcfc12    )
         deallocate (colcfc22    )
         deallocate (colccl4     )
         deallocate (plankbnd    )
         deallocate (dplankbnd_dt)
         deallocate (planklev    )
         deallocate (planklay    )
         deallocate (jp          )
         deallocate (jt          )
         deallocate (jt1         )
         deallocate (forfac      )
         deallocate (indfor      )
         deallocate (forfrac     )
         deallocate (selffac     )
         deallocate (indself     )
         deallocate (selffrac    )
         deallocate (scaleminor  )
         deallocate (scaleminorn2)
         deallocate (indminor    )
         deallocate (minorfrac   )
         deallocate (rat_h2oco2  )
         deallocate (rat_h2oco2_1)
         deallocate (rat_h2oo3   )
         deallocate (rat_h2oo3_1 )
         deallocate (rat_h2on2o  )
         deallocate (rat_h2on2o_1)
         deallocate (rat_h2och4  )
         deallocate (rat_h2och4_1)
         deallocate (rat_n2oco2  )
         deallocate (rat_n2oco2_1)
         deallocate (rat_o3co2   )
         deallocate (rat_o3co2_1 )
         deallocate (fac00       )
         deallocate (fac01       )
         deallocate (fac10       )
         deallocate (fac11       )
      end if
   end subroutine setcoef_free


!***************************************************************************
   subroutine lwatmref
!***************************************************************************

      save
 
! These pressures are chosen such that the ln of the first pressure
! has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
! each subsequent ln(pressure) differs from the previous one by 0.2.

      pref(:) = (/ &
          1.05363e+03 ,8.62642e+02 ,7.06272e+02 ,5.78246e+02 ,4.73428e+02 , &
          3.87610e+02 ,3.17348e+02 ,2.59823e+02 ,2.12725e+02 ,1.74164e+02 , &
          1.42594e+02 ,1.16746e+02 ,9.55835e+01 ,7.82571e+01 ,6.40715e+01 , &
          5.24573e+01 ,4.29484e+01 ,3.51632e+01 ,2.87892e+01 ,2.35706e+01 , &
          1.92980e+01 ,1.57998e+01 ,1.29358e+01 ,1.05910e+01 ,8.67114e+00 , &
          7.09933e+00 ,5.81244e+00 ,4.75882e+00 ,3.89619e+00 ,3.18993e+00 , &
          2.61170e+00 ,2.13828e+00 ,1.75067e+00 ,1.43333e+00 ,1.17351e+00 , &
          9.60789e-01 ,7.86628e-01 ,6.44036e-01 ,5.27292e-01 ,4.31710e-01 , &
          3.53455e-01 ,2.89384e-01 ,2.36928e-01 ,1.93980e-01 ,1.58817e-01 , &
          1.30029e-01 ,1.06458e-01 ,8.71608e-02 ,7.13612e-02 ,5.84256e-02 , &
          4.78349e-02 ,3.91639e-02 ,3.20647e-02 ,2.62523e-02 ,2.14936e-02 , &
          1.75975e-02 ,1.44076e-02 ,1.17959e-02 ,9.65769e-03 /)

      preflog(:) = (/ &
           6.9600e+00 , 6.7600e+00 , 6.5600e+00 , 6.3600e+00 , 6.1600e+00 , &
           5.9600e+00 , 5.7600e+00 , 5.5600e+00 , 5.3600e+00 , 5.1600e+00 , &
           4.9600e+00 , 4.7600e+00 , 4.5600e+00 , 4.3600e+00 , 4.1600e+00 , &
           3.9600e+00 , 3.7600e+00 , 3.5600e+00 , 3.3600e+00 , 3.1600e+00 , &
           2.9600e+00 , 2.7600e+00 , 2.5600e+00 , 2.3600e+00 , 2.1600e+00 , &
           1.9600e+00 , 1.7600e+00 , 1.5600e+00 , 1.3600e+00 , 1.1600e+00 , &
           9.6000e-01 , 7.6000e-01 , 5.6000e-01 , 3.6000e-01 , 1.6000e-01 , &
          -4.0000e-02 ,-2.4000e-01 ,-4.4000e-01 ,-6.4000e-01 ,-8.4000e-01 , &
          -1.0400e+00 ,-1.2400e+00 ,-1.4400e+00 ,-1.6400e+00 ,-1.8400e+00 , &
          -2.0400e+00 ,-2.2400e+00 ,-2.4400e+00 ,-2.6400e+00 ,-2.8400e+00 , &
          -3.0400e+00 ,-3.2400e+00 ,-3.4400e+00 ,-3.6400e+00 ,-3.8400e+00 , &
          -4.0400e+00 ,-4.2400e+00 ,-4.4400e+00 ,-4.6400e+00 /)

! These are the temperatures associated with the respective 
! pressures for the mls standard atmosphere. 

      tref(:) = (/ &
           2.9420e+02 , 2.8799e+02 , 2.7894e+02 , 2.6925e+02 , 2.5983e+02 , &
           2.5017e+02 , 2.4077e+02 , 2.3179e+02 , 2.2306e+02 , 2.1578e+02 , &
           2.1570e+02 , 2.1570e+02 , 2.1570e+02 , 2.1706e+02 , 2.1858e+02 , &
           2.2018e+02 , 2.2174e+02 , 2.2328e+02 , 2.2479e+02 , 2.2655e+02 , &
           2.2834e+02 , 2.3113e+02 , 2.3401e+02 , 2.3703e+02 , 2.4022e+02 , &
           2.4371e+02 , 2.4726e+02 , 2.5085e+02 , 2.5457e+02 , 2.5832e+02 , &
           2.6216e+02 , 2.6606e+02 , 2.6999e+02 , 2.7340e+02 , 2.7536e+02 , &
           2.7568e+02 , 2.7372e+02 , 2.7163e+02 , 2.6955e+02 , 2.6593e+02 , &
           2.6211e+02 , 2.5828e+02 , 2.5360e+02 , 2.4854e+02 , 2.4348e+02 , &
           2.3809e+02 , 2.3206e+02 , 2.2603e+02 , 2.2000e+02 , 2.1435e+02 , &
           2.0887e+02 , 2.0340e+02 , 1.9792e+02 , 1.9290e+02 , 1.8809e+02 , &
           1.8329e+02 , 1.7849e+02 , 1.7394e+02 , 1.7212e+02 /)

       chi_mls(1,1:12) = (/ &
        1.8760e-02 , 1.2223e-02 , 5.8909e-03 , 2.7675e-03 , 1.4065e-03 , &
        7.5970e-04 , 3.8876e-04 , 1.6542e-04 , 3.7190e-05 , 7.4765e-06 , &
        4.3082e-06 , 3.3319e-06 /)
       chi_mls(1,13:59) = (/ &
        3.2039e-06 ,  3.1619e-06 ,  3.2524e-06 ,  3.4226e-06 ,  3.6288e-06 , &
        3.9148e-06 ,  4.1488e-06 ,  4.3081e-06 ,  4.4420e-06 ,  4.5778e-06 , &
        4.7087e-06 ,  4.7943e-06 ,  4.8697e-06 ,  4.9260e-06 ,  4.9669e-06 , &
        4.9963e-06 ,  5.0527e-06 ,  5.1266e-06 ,  5.2503e-06 ,  5.3571e-06 , &
        5.4509e-06 ,  5.4830e-06 ,  5.5000e-06 ,  5.5000e-06 ,  5.4536e-06 , &
        5.4047e-06 ,  5.3558e-06 ,  5.2533e-06 ,  5.1436e-06 ,  5.0340e-06 , &
        4.8766e-06 ,  4.6979e-06 ,  4.5191e-06 ,  4.3360e-06 ,  4.1442e-06 , &
        3.9523e-06 ,  3.7605e-06 ,  3.5722e-06 ,  3.3855e-06 ,  3.1988e-06 , &
        3.0121e-06 ,  2.8262e-06 ,  2.6407e-06 ,  2.4552e-06 ,  2.2696e-06 , &
        4.3360e-06 ,  4.1442e-06 /)
       chi_mls(2,1:12) = (/ &
        3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 , &
        3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 , &
        3.5500e-04 ,  3.5500e-04 /)
       chi_mls(2,13:59) = (/ &
        3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 , &
        3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 , &
        3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 , &
        3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 , &
        3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 , &
        3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 , &
        3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 , &
        3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 ,  3.5500e-04 , &
        3.5500e-04 ,  3.5471e-04 ,  3.5427e-04 ,  3.5384e-04 ,  3.5340e-04 , &
        3.5500e-04 ,  3.5500e-04 /)
       chi_mls(3,1:12) = (/ &
        3.0170e-08 ,  3.4725e-08 ,  4.2477e-08 ,  5.2759e-08 ,  6.6944e-08 , &
        8.7130e-08 ,  1.1391e-07 ,  1.5677e-07 ,  2.1788e-07 ,  3.2443e-07 , &
        4.6594e-07 ,  5.6806e-07 /)
       chi_mls(3,13:59) = (/ &
        6.9607e-07 ,  1.1186e-06 ,  1.7618e-06 ,  2.3269e-06 ,  2.9577e-06 , &
        3.6593e-06 ,  4.5950e-06 ,  5.3189e-06 ,  5.9618e-06 ,  6.5113e-06 , &
        7.0635e-06 ,  7.6917e-06 ,  8.2577e-06 ,  8.7082e-06 ,  8.8325e-06 , &
        8.7149e-06 ,  8.0943e-06 ,  7.3307e-06 ,  6.3101e-06 ,  5.3672e-06 , &
        4.4829e-06 ,  3.8391e-06 ,  3.2827e-06 ,  2.8235e-06 ,  2.4906e-06 , &
        2.1645e-06 ,  1.8385e-06 ,  1.6618e-06 ,  1.5052e-06 ,  1.3485e-06 , &
        1.1972e-06 ,  1.0482e-06 ,  8.9926e-07 ,  7.6343e-07 ,  6.5381e-07 , &
        5.4419e-07 ,  4.3456e-07 ,  3.6421e-07 ,  3.1194e-07 ,  2.5967e-07 , &
        2.0740e-07 ,  1.9146e-07 ,  1.9364e-07 ,  1.9582e-07 ,  1.9800e-07 , &
        7.6343e-07 ,  6.5381e-07 /)
       chi_mls(4,1:12) = (/ &
        3.2000e-07 ,  3.2000e-07 ,  3.2000e-07 ,  3.2000e-07 ,  3.2000e-07 , &
        3.1965e-07 ,  3.1532e-07 ,  3.0383e-07 ,  2.9422e-07 ,  2.8495e-07 , &
        2.7671e-07 ,  2.6471e-07 /)
       chi_mls(4,13:59) = (/ &
        2.4285e-07 ,  2.0955e-07 ,  1.7195e-07 ,  1.3749e-07 ,  1.1332e-07 , &
        1.0035e-07 ,  9.1281e-08 ,  8.5463e-08 ,  8.0363e-08 ,  7.3372e-08 , &
        6.5975e-08 ,  5.6039e-08 ,  4.7090e-08 ,  3.9977e-08 ,  3.2979e-08 , &
        2.6064e-08 ,  2.1066e-08 ,  1.6592e-08 ,  1.3017e-08 ,  1.0090e-08 , &
        7.6249e-09 ,  6.1159e-09 ,  4.6672e-09 ,  3.2857e-09 ,  2.8484e-09 , &
        2.4620e-09 ,  2.0756e-09 ,  1.8551e-09 ,  1.6568e-09 ,  1.4584e-09 , &
        1.3195e-09 ,  1.2072e-09 ,  1.0948e-09 ,  9.9780e-10 ,  9.3126e-10 , &
        8.6472e-10 ,  7.9818e-10 ,  7.5138e-10 ,  7.1367e-10 ,  6.7596e-10 , &
        6.3825e-10 ,  6.0981e-10 ,  5.8600e-10 ,  5.6218e-10 ,  5.3837e-10 , &
        9.9780e-10 ,  9.3126e-10 /)
       chi_mls(5,1:12) = (/ &
        1.5000e-07 ,  1.4306e-07 ,  1.3474e-07 ,  1.3061e-07 ,  1.2793e-07 , &
        1.2038e-07 ,  1.0798e-07 ,  9.4238e-08 ,  7.9488e-08 ,  6.1386e-08 , &
        4.5563e-08 ,  3.3475e-08 /)
       chi_mls(5,13:59) = (/ &
        2.5118e-08 ,  1.8671e-08 ,  1.4349e-08 ,  1.2501e-08 ,  1.2407e-08 , &
        1.3472e-08 ,  1.4900e-08 ,  1.6079e-08 ,  1.7156e-08 ,  1.8616e-08 , &
        2.0106e-08 ,  2.1654e-08 ,  2.3096e-08 ,  2.4340e-08 ,  2.5643e-08 , &
        2.6990e-08 ,  2.8456e-08 ,  2.9854e-08 ,  3.0943e-08 ,  3.2023e-08 , &
        3.3101e-08 ,  3.4260e-08 ,  3.5360e-08 ,  3.6397e-08 ,  3.7310e-08 , &
        3.8217e-08 ,  3.9123e-08 ,  4.1303e-08 ,  4.3652e-08 ,  4.6002e-08 , &
        5.0289e-08 ,  5.5446e-08 ,  6.0603e-08 ,  6.8946e-08 ,  8.3652e-08 , &
        9.8357e-08 ,  1.1306e-07 ,  1.4766e-07 ,  1.9142e-07 ,  2.3518e-07 , &
        2.7894e-07 ,  3.5001e-07 ,  4.3469e-07 ,  5.1938e-07 ,  6.0407e-07 , &
        6.8946e-08 ,  8.3652e-08 /)
       chi_mls(6,1:12) = (/ &
        1.7000e-06 ,  1.7000e-06 ,  1.6999e-06 ,  1.6904e-06 ,  1.6671e-06 , &
        1.6351e-06 ,  1.6098e-06 ,  1.5590e-06 ,  1.5120e-06 ,  1.4741e-06 , &
        1.4385e-06 ,  1.4002e-06 /)
       chi_mls(6,13:59) = (/ &
        1.3573e-06 ,  1.3130e-06 ,  1.2512e-06 ,  1.1668e-06 ,  1.0553e-06 , &
        9.3281e-07 ,  8.1217e-07 ,  7.5239e-07 ,  7.0728e-07 ,  6.6722e-07 , &
        6.2733e-07 ,  5.8604e-07 ,  5.4769e-07 ,  5.1480e-07 ,  4.8206e-07 , &
        4.4943e-07 ,  4.1702e-07 ,  3.8460e-07 ,  3.5200e-07 ,  3.1926e-07 , &
        2.8646e-07 ,  2.5498e-07 ,  2.2474e-07 ,  1.9588e-07 ,  1.8295e-07 , &
        1.7089e-07 ,  1.5882e-07 ,  1.5536e-07 ,  1.5304e-07 ,  1.5072e-07 , &
        1.5000e-07 ,  1.5000e-07 ,  1.5000e-07 ,  1.5000e-07 ,  1.5000e-07 , &
        1.5000e-07 ,  1.5000e-07 ,  1.5000e-07 ,  1.5000e-07 ,  1.5000e-07 , &
        1.5000e-07 ,  1.5000e-07 ,  1.5000e-07 ,  1.5000e-07 ,  1.5000e-07 , &
        1.5000e-07 ,  1.5000e-07 /)
       chi_mls(7,1:12) = (/ &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 /)
       chi_mls(7,13:59) = (/ &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 ,  0.2090 ,  0.2090 ,  0.2090 , &
        0.2090 ,  0.2090 /)

   end subroutine lwatmref

!***************************************************************************
   subroutine lwavplank
!***************************************************************************

      save
 
      totplnk(1:50,  1) = (/ &
      0.14783e-05 ,0.15006e-05 ,0.15230e-05 ,0.15455e-05 ,0.15681e-05 , &
      0.15908e-05 ,0.16136e-05 ,0.16365e-05 ,0.16595e-05 ,0.16826e-05 , &
      0.17059e-05 ,0.17292e-05 ,0.17526e-05 ,0.17762e-05 ,0.17998e-05 , &
      0.18235e-05 ,0.18473e-05 ,0.18712e-05 ,0.18953e-05 ,0.19194e-05 , &
      0.19435e-05 ,0.19678e-05 ,0.19922e-05 ,0.20166e-05 ,0.20412e-05 , &
      0.20658e-05 ,0.20905e-05 ,0.21153e-05 ,0.21402e-05 ,0.21652e-05 , &
      0.21902e-05 ,0.22154e-05 ,0.22406e-05 ,0.22659e-05 ,0.22912e-05 , &
      0.23167e-05 ,0.23422e-05 ,0.23678e-05 ,0.23934e-05 ,0.24192e-05 , &
      0.24450e-05 ,0.24709e-05 ,0.24968e-05 ,0.25229e-05 ,0.25490e-05 , &
      0.25751e-05 ,0.26014e-05 ,0.26277e-05 ,0.26540e-05 ,0.26805e-05 /)
      totplnk(51:100,  1) = (/ &
      0.27070e-05 ,0.27335e-05 ,0.27602e-05 ,0.27869e-05 ,0.28136e-05 , &
      0.28404e-05 ,0.28673e-05 ,0.28943e-05 ,0.29213e-05 ,0.29483e-05 , &
      0.29754e-05 ,0.30026e-05 ,0.30298e-05 ,0.30571e-05 ,0.30845e-05 , &
      0.31119e-05 ,0.31393e-05 ,0.31669e-05 ,0.31944e-05 ,0.32220e-05 , &
      0.32497e-05 ,0.32774e-05 ,0.33052e-05 ,0.33330e-05 ,0.33609e-05 , &
      0.33888e-05 ,0.34168e-05 ,0.34448e-05 ,0.34729e-05 ,0.35010e-05 , &
      0.35292e-05 ,0.35574e-05 ,0.35857e-05 ,0.36140e-05 ,0.36424e-05 , &
      0.36708e-05 ,0.36992e-05 ,0.37277e-05 ,0.37563e-05 ,0.37848e-05 , &
      0.38135e-05 ,0.38421e-05 ,0.38708e-05 ,0.38996e-05 ,0.39284e-05 , &
      0.39572e-05 ,0.39861e-05 ,0.40150e-05 ,0.40440e-05 ,0.40730e-05 /)
      totplnk(101:150,  1) = (/ &
      0.41020e-05 ,0.41311e-05 ,0.41602e-05 ,0.41893e-05 ,0.42185e-05 , &
      0.42477e-05 ,0.42770e-05 ,0.43063e-05 ,0.43356e-05 ,0.43650e-05 , &
      0.43944e-05 ,0.44238e-05 ,0.44533e-05 ,0.44828e-05 ,0.45124e-05 , &
      0.45419e-05 ,0.45715e-05 ,0.46012e-05 ,0.46309e-05 ,0.46606e-05 , &
      0.46903e-05 ,0.47201e-05 ,0.47499e-05 ,0.47797e-05 ,0.48096e-05 , &
      0.48395e-05 ,0.48695e-05 ,0.48994e-05 ,0.49294e-05 ,0.49594e-05 , &
      0.49895e-05 ,0.50196e-05 ,0.50497e-05 ,0.50798e-05 ,0.51100e-05 , &
      0.51402e-05 ,0.51704e-05 ,0.52007e-05 ,0.52309e-05 ,0.52612e-05 , &
      0.52916e-05 ,0.53219e-05 ,0.53523e-05 ,0.53827e-05 ,0.54132e-05 , &
      0.54436e-05 ,0.54741e-05 ,0.55047e-05 ,0.55352e-05 ,0.55658e-05 /)
      totplnk(151:181,  1) = (/ &
      0.55964e-05 ,0.56270e-05 ,0.56576e-05 ,0.56883e-05 ,0.57190e-05 , &
      0.57497e-05 ,0.57804e-05 ,0.58112e-05 ,0.58420e-05 ,0.58728e-05 , &
      0.59036e-05 ,0.59345e-05 ,0.59653e-05 ,0.59962e-05 ,0.60272e-05 , &
      0.60581e-05 ,0.60891e-05 ,0.61201e-05 ,0.61511e-05 ,0.61821e-05 , &
      0.62131e-05 ,0.62442e-05 ,0.62753e-05 ,0.63064e-05 ,0.63376e-05 , &
      0.63687e-05 ,0.63998e-05 ,0.64310e-05 ,0.64622e-05 ,0.64935e-05 , &
      0.65247e-05 /)
      totplnk(1:50,  2) = (/ &
      0.20262e-05 ,0.20757e-05 ,0.21257e-05 ,0.21763e-05 ,0.22276e-05 , &
      0.22794e-05 ,0.23319e-05 ,0.23849e-05 ,0.24386e-05 ,0.24928e-05 , &
      0.25477e-05 ,0.26031e-05 ,0.26591e-05 ,0.27157e-05 ,0.27728e-05 , &
      0.28306e-05 ,0.28889e-05 ,0.29478e-05 ,0.30073e-05 ,0.30673e-05 , &
      0.31279e-05 ,0.31890e-05 ,0.32507e-05 ,0.33129e-05 ,0.33757e-05 , &
      0.34391e-05 ,0.35029e-05 ,0.35674e-05 ,0.36323e-05 ,0.36978e-05 , &
      0.37638e-05 ,0.38304e-05 ,0.38974e-05 ,0.39650e-05 ,0.40331e-05 , &
      0.41017e-05 ,0.41708e-05 ,0.42405e-05 ,0.43106e-05 ,0.43812e-05 , &
      0.44524e-05 ,0.45240e-05 ,0.45961e-05 ,0.46687e-05 ,0.47418e-05 , &
      0.48153e-05 ,0.48894e-05 ,0.49639e-05 ,0.50389e-05 ,0.51143e-05 /)
      totplnk(51:100,  2) = (/ &
      0.51902e-05 ,0.52666e-05 ,0.53434e-05 ,0.54207e-05 ,0.54985e-05 , &
      0.55767e-05 ,0.56553e-05 ,0.57343e-05 ,0.58139e-05 ,0.58938e-05 , &
      0.59742e-05 ,0.60550e-05 ,0.61362e-05 ,0.62179e-05 ,0.63000e-05 , &
      0.63825e-05 ,0.64654e-05 ,0.65487e-05 ,0.66324e-05 ,0.67166e-05 , &
      0.68011e-05 ,0.68860e-05 ,0.69714e-05 ,0.70571e-05 ,0.71432e-05 , &
      0.72297e-05 ,0.73166e-05 ,0.74039e-05 ,0.74915e-05 ,0.75796e-05 , &
      0.76680e-05 ,0.77567e-05 ,0.78459e-05 ,0.79354e-05 ,0.80252e-05 , &
      0.81155e-05 ,0.82061e-05 ,0.82970e-05 ,0.83883e-05 ,0.84799e-05 , &
      0.85719e-05 ,0.86643e-05 ,0.87569e-05 ,0.88499e-05 ,0.89433e-05 , &
      0.90370e-05 ,0.91310e-05 ,0.92254e-05 ,0.93200e-05 ,0.94150e-05 /)
      totplnk(101:150,  2) = (/ &
      0.95104e-05 ,0.96060e-05 ,0.97020e-05 ,0.97982e-05 ,0.98948e-05 , &
      0.99917e-05 ,0.10089e-04 ,0.10186e-04 ,0.10284e-04 ,0.10382e-04 , &
      0.10481e-04 ,0.10580e-04 ,0.10679e-04 ,0.10778e-04 ,0.10877e-04 , &
      0.10977e-04 ,0.11077e-04 ,0.11178e-04 ,0.11279e-04 ,0.11380e-04 , &
      0.11481e-04 ,0.11583e-04 ,0.11684e-04 ,0.11786e-04 ,0.11889e-04 , &
      0.11992e-04 ,0.12094e-04 ,0.12198e-04 ,0.12301e-04 ,0.12405e-04 , &
      0.12509e-04 ,0.12613e-04 ,0.12717e-04 ,0.12822e-04 ,0.12927e-04 , &
      0.13032e-04 ,0.13138e-04 ,0.13244e-04 ,0.13349e-04 ,0.13456e-04 , &
      0.13562e-04 ,0.13669e-04 ,0.13776e-04 ,0.13883e-04 ,0.13990e-04 , &
      0.14098e-04 ,0.14206e-04 ,0.14314e-04 ,0.14422e-04 ,0.14531e-04 /)
      totplnk(151:181,  2) = (/ &
      0.14639e-04 ,0.14748e-04 ,0.14857e-04 ,0.14967e-04 ,0.15076e-04 , &
      0.15186e-04 ,0.15296e-04 ,0.15407e-04 ,0.15517e-04 ,0.15628e-04 , &
      0.15739e-04 ,0.15850e-04 ,0.15961e-04 ,0.16072e-04 ,0.16184e-04 , &
      0.16296e-04 ,0.16408e-04 ,0.16521e-04 ,0.16633e-04 ,0.16746e-04 , &
      0.16859e-04 ,0.16972e-04 ,0.17085e-04 ,0.17198e-04 ,0.17312e-04 , &
      0.17426e-04 ,0.17540e-04 ,0.17654e-04 ,0.17769e-04 ,0.17883e-04 , &
      0.17998e-04 /)
      totplnk(1:50, 3) = (/ &
      1.34822e-06 ,1.39134e-06 ,1.43530e-06 ,1.48010e-06 ,1.52574e-06 , &
      1.57222e-06 ,1.61956e-06 ,1.66774e-06 ,1.71678e-06 ,1.76666e-06 , &
      1.81741e-06 ,1.86901e-06 ,1.92147e-06 ,1.97479e-06 ,2.02898e-06 , &
      2.08402e-06 ,2.13993e-06 ,2.19671e-06 ,2.25435e-06 ,2.31285e-06 , &
      2.37222e-06 ,2.43246e-06 ,2.49356e-06 ,2.55553e-06 ,2.61837e-06 , &
      2.68207e-06 ,2.74664e-06 ,2.81207e-06 ,2.87837e-06 ,2.94554e-06 , &
      3.01356e-06 ,3.08245e-06 ,3.15221e-06 ,3.22282e-06 ,3.29429e-06 , &
      3.36662e-06 ,3.43982e-06 ,3.51386e-06 ,3.58876e-06 ,3.66451e-06 , &
      3.74112e-06 ,3.81857e-06 ,3.89688e-06 ,3.97602e-06 ,4.05601e-06 , &
      4.13685e-06 ,4.21852e-06 ,4.30104e-06 ,4.38438e-06 ,4.46857e-06 /)
      totplnk(51:100, 3) = (/ &
      4.55358e-06 ,4.63943e-06 ,4.72610e-06 ,4.81359e-06 ,4.90191e-06 , &
      4.99105e-06 ,5.08100e-06 ,5.17176e-06 ,5.26335e-06 ,5.35573e-06 , &
      5.44892e-06 ,5.54292e-06 ,5.63772e-06 ,5.73331e-06 ,5.82970e-06 , &
      5.92688e-06 ,6.02485e-06 ,6.12360e-06 ,6.22314e-06 ,6.32346e-06 , &
      6.42455e-06 ,6.52641e-06 ,6.62906e-06 ,6.73247e-06 ,6.83664e-06 , &
      6.94156e-06 ,7.04725e-06 ,7.15370e-06 ,7.26089e-06 ,7.36883e-06 , &
      7.47752e-06 ,7.58695e-06 ,7.69712e-06 ,7.80801e-06 ,7.91965e-06 , &
      8.03201e-06 ,8.14510e-06 ,8.25891e-06 ,8.37343e-06 ,8.48867e-06 , &
      8.60463e-06 ,8.72128e-06 ,8.83865e-06 ,8.95672e-06 ,9.07548e-06 , &
      9.19495e-06 ,9.31510e-06 ,9.43594e-06 ,9.55745e-06 ,9.67966e-06 /)
      totplnk(101:150, 3) = (/ &
      9.80254e-06 ,9.92609e-06 ,1.00503e-05 ,1.01752e-05 ,1.03008e-05 , &
      1.04270e-05 ,1.05539e-05 ,1.06814e-05 ,1.08096e-05 ,1.09384e-05 , &
      1.10679e-05 ,1.11980e-05 ,1.13288e-05 ,1.14601e-05 ,1.15922e-05 , &
      1.17248e-05 ,1.18581e-05 ,1.19920e-05 ,1.21265e-05 ,1.22616e-05 , &
      1.23973e-05 ,1.25337e-05 ,1.26706e-05 ,1.28081e-05 ,1.29463e-05 , &
      1.30850e-05 ,1.32243e-05 ,1.33642e-05 ,1.35047e-05 ,1.36458e-05 , &
      1.37875e-05 ,1.39297e-05 ,1.40725e-05 ,1.42159e-05 ,1.43598e-05 , &
      1.45044e-05 ,1.46494e-05 ,1.47950e-05 ,1.49412e-05 ,1.50879e-05 , &
      1.52352e-05 ,1.53830e-05 ,1.55314e-05 ,1.56803e-05 ,1.58297e-05 , &
      1.59797e-05 ,1.61302e-05 ,1.62812e-05 ,1.64327e-05 ,1.65848e-05 /)
      totplnk(151:181, 3) = (/ &
      1.67374e-05 ,1.68904e-05 ,1.70441e-05 ,1.71982e-05 ,1.73528e-05 , &
      1.75079e-05 ,1.76635e-05 ,1.78197e-05 ,1.79763e-05 ,1.81334e-05 , &
      1.82910e-05 ,1.84491e-05 ,1.86076e-05 ,1.87667e-05 ,1.89262e-05 , &
      1.90862e-05 ,1.92467e-05 ,1.94076e-05 ,1.95690e-05 ,1.97309e-05 , &
      1.98932e-05 ,2.00560e-05 ,2.02193e-05 ,2.03830e-05 ,2.05472e-05 , &
      2.07118e-05 ,2.08768e-05 ,2.10423e-05 ,2.12083e-05 ,2.13747e-05 , &
      2.15414e-05 /)
      totplnk(1:50, 4) = (/ &
      8.90528e-07 ,9.24222e-07 ,9.58757e-07 ,9.94141e-07 ,1.03038e-06 , &
      1.06748e-06 ,1.10545e-06 ,1.14430e-06 ,1.18403e-06 ,1.22465e-06 , &
      1.26618e-06 ,1.30860e-06 ,1.35193e-06 ,1.39619e-06 ,1.44136e-06 , &
      1.48746e-06 ,1.53449e-06 ,1.58246e-06 ,1.63138e-06 ,1.68124e-06 , &
      1.73206e-06 ,1.78383e-06 ,1.83657e-06 ,1.89028e-06 ,1.94495e-06 , &
      2.00060e-06 ,2.05724e-06 ,2.11485e-06 ,2.17344e-06 ,2.23303e-06 , &
      2.29361e-06 ,2.35519e-06 ,2.41777e-06 ,2.48134e-06 ,2.54592e-06 , &
      2.61151e-06 ,2.67810e-06 ,2.74571e-06 ,2.81433e-06 ,2.88396e-06 , &
      2.95461e-06 ,3.02628e-06 ,3.09896e-06 ,3.17267e-06 ,3.24741e-06 , &
      3.32316e-06 ,3.39994e-06 ,3.47774e-06 ,3.55657e-06 ,3.63642e-06 /)
      totplnk(51:100, 4) = (/ &
      3.71731e-06 ,3.79922e-06 ,3.88216e-06 ,3.96612e-06 ,4.05112e-06 , &
      4.13714e-06 ,4.22419e-06 ,4.31227e-06 ,4.40137e-06 ,4.49151e-06 , &
      4.58266e-06 ,4.67485e-06 ,4.76806e-06 ,4.86229e-06 ,4.95754e-06 , &
      5.05383e-06 ,5.15113e-06 ,5.24946e-06 ,5.34879e-06 ,5.44916e-06 , &
      5.55053e-06 ,5.65292e-06 ,5.75632e-06 ,5.86073e-06 ,5.96616e-06 , &
      6.07260e-06 ,6.18003e-06 ,6.28848e-06 ,6.39794e-06 ,6.50838e-06 , &
      6.61983e-06 ,6.73229e-06 ,6.84573e-06 ,6.96016e-06 ,7.07559e-06 , &
      7.19200e-06 ,7.30940e-06 ,7.42779e-06 ,7.54715e-06 ,7.66749e-06 , &
      7.78882e-06 ,7.91110e-06 ,8.03436e-06 ,8.15859e-06 ,8.28379e-06 , &
      8.40994e-06 ,8.53706e-06 ,8.66515e-06 ,8.79418e-06 ,8.92416e-06 /)
      totplnk(101:150, 4) = (/ &
      9.05510e-06 ,9.18697e-06 ,9.31979e-06 ,9.45356e-06 ,9.58826e-06 , &
      9.72389e-06 ,9.86046e-06 ,9.99793e-06 ,1.01364e-05 ,1.02757e-05 , &
      1.04159e-05 ,1.05571e-05 ,1.06992e-05 ,1.08422e-05 ,1.09861e-05 , &
      1.11309e-05 ,1.12766e-05 ,1.14232e-05 ,1.15707e-05 ,1.17190e-05 , &
      1.18683e-05 ,1.20184e-05 ,1.21695e-05 ,1.23214e-05 ,1.24741e-05 , &
      1.26277e-05 ,1.27822e-05 ,1.29376e-05 ,1.30939e-05 ,1.32509e-05 , &
      1.34088e-05 ,1.35676e-05 ,1.37273e-05 ,1.38877e-05 ,1.40490e-05 , &
      1.42112e-05 ,1.43742e-05 ,1.45380e-05 ,1.47026e-05 ,1.48680e-05 , &
      1.50343e-05 ,1.52014e-05 ,1.53692e-05 ,1.55379e-05 ,1.57074e-05 , &
      1.58778e-05 ,1.60488e-05 ,1.62207e-05 ,1.63934e-05 ,1.65669e-05 /)
      totplnk(151:181, 4) = (/ &
      1.67411e-05 ,1.69162e-05 ,1.70920e-05 ,1.72685e-05 ,1.74459e-05 , &
      1.76240e-05 ,1.78029e-05 ,1.79825e-05 ,1.81629e-05 ,1.83440e-05 , &
      1.85259e-05 ,1.87086e-05 ,1.88919e-05 ,1.90760e-05 ,1.92609e-05 , &
      1.94465e-05 ,1.96327e-05 ,1.98199e-05 ,2.00076e-05 ,2.01961e-05 , &
      2.03853e-05 ,2.05752e-05 ,2.07658e-05 ,2.09571e-05 ,2.11491e-05 , &
      2.13418e-05 ,2.15352e-05 ,2.17294e-05 ,2.19241e-05 ,2.21196e-05 , &
      2.23158e-05 /)
      totplnk(1:50, 5) = (/ &
      5.70230e-07 ,5.94788e-07 ,6.20085e-07 ,6.46130e-07 ,6.72936e-07 , &
      7.00512e-07 ,7.28869e-07 ,7.58019e-07 ,7.87971e-07 ,8.18734e-07 , &
      8.50320e-07 ,8.82738e-07 ,9.15999e-07 ,9.50110e-07 ,9.85084e-07 , &
      1.02093e-06 ,1.05765e-06 ,1.09527e-06 ,1.13378e-06 ,1.17320e-06 , &
      1.21353e-06 ,1.25479e-06 ,1.29698e-06 ,1.34011e-06 ,1.38419e-06 , &
      1.42923e-06 ,1.47523e-06 ,1.52221e-06 ,1.57016e-06 ,1.61910e-06 , &
      1.66904e-06 ,1.71997e-06 ,1.77192e-06 ,1.82488e-06 ,1.87886e-06 , &
      1.93387e-06 ,1.98991e-06 ,2.04699e-06 ,2.10512e-06 ,2.16430e-06 , &
      2.22454e-06 ,2.28584e-06 ,2.34821e-06 ,2.41166e-06 ,2.47618e-06 , &
      2.54178e-06 ,2.60847e-06 ,2.67626e-06 ,2.74514e-06 ,2.81512e-06 /)
      totplnk(51:100, 5) = (/ &
      2.88621e-06 ,2.95841e-06 ,3.03172e-06 ,3.10615e-06 ,3.18170e-06 , &
      3.25838e-06 ,3.33618e-06 ,3.41511e-06 ,3.49518e-06 ,3.57639e-06 , &
      3.65873e-06 ,3.74221e-06 ,3.82684e-06 ,3.91262e-06 ,3.99955e-06 , &
      4.08763e-06 ,4.17686e-06 ,4.26725e-06 ,4.35880e-06 ,4.45150e-06 , &
      4.54537e-06 ,4.64039e-06 ,4.73659e-06 ,4.83394e-06 ,4.93246e-06 , &
      5.03215e-06 ,5.13301e-06 ,5.23504e-06 ,5.33823e-06 ,5.44260e-06 , &
      5.54814e-06 ,5.65484e-06 ,5.76272e-06 ,5.87177e-06 ,5.98199e-06 , &
      6.09339e-06 ,6.20596e-06 ,6.31969e-06 ,6.43460e-06 ,6.55068e-06 , &
      6.66793e-06 ,6.78636e-06 ,6.90595e-06 ,7.02670e-06 ,7.14863e-06 , &
      7.27173e-06 ,7.39599e-06 ,7.52142e-06 ,7.64802e-06 ,7.77577e-06 /)
      totplnk(101:150, 5) = (/ &
      7.90469e-06 ,8.03477e-06 ,8.16601e-06 ,8.29841e-06 ,8.43198e-06 , &
      8.56669e-06 ,8.70256e-06 ,8.83957e-06 ,8.97775e-06 ,9.11706e-06 , &
      9.25753e-06 ,9.39915e-06 ,9.54190e-06 ,9.68580e-06 ,9.83085e-06 , &
      9.97704e-06 ,1.01243e-05 ,1.02728e-05 ,1.04224e-05 ,1.05731e-05 , &
      1.07249e-05 ,1.08779e-05 ,1.10320e-05 ,1.11872e-05 ,1.13435e-05 , &
      1.15009e-05 ,1.16595e-05 ,1.18191e-05 ,1.19799e-05 ,1.21418e-05 , &
      1.23048e-05 ,1.24688e-05 ,1.26340e-05 ,1.28003e-05 ,1.29676e-05 , &
      1.31361e-05 ,1.33056e-05 ,1.34762e-05 ,1.36479e-05 ,1.38207e-05 , &
      1.39945e-05 ,1.41694e-05 ,1.43454e-05 ,1.45225e-05 ,1.47006e-05 , &
      1.48797e-05 ,1.50600e-05 ,1.52413e-05 ,1.54236e-05 ,1.56070e-05 /)
      totplnk(151:181, 5) = (/ &
      1.57914e-05 ,1.59768e-05 ,1.61633e-05 ,1.63509e-05 ,1.65394e-05 , &
      1.67290e-05 ,1.69197e-05 ,1.71113e-05 ,1.73040e-05 ,1.74976e-05 , &
      1.76923e-05 ,1.78880e-05 ,1.80847e-05 ,1.82824e-05 ,1.84811e-05 , &
      1.86808e-05 ,1.88814e-05 ,1.90831e-05 ,1.92857e-05 ,1.94894e-05 , &
      1.96940e-05 ,1.98996e-05 ,2.01061e-05 ,2.03136e-05 ,2.05221e-05 , &
      2.07316e-05 ,2.09420e-05 ,2.11533e-05 ,2.13657e-05 ,2.15789e-05 , &
      2.17931e-05 /)
      totplnk(1:50, 6) = (/ &
      2.73493e-07 ,2.87408e-07 ,3.01848e-07 ,3.16825e-07 ,3.32352e-07 , &
      3.48439e-07 ,3.65100e-07 ,3.82346e-07 ,4.00189e-07 ,4.18641e-07 , &
      4.37715e-07 ,4.57422e-07 ,4.77774e-07 ,4.98784e-07 ,5.20464e-07 , &
      5.42824e-07 ,5.65879e-07 ,5.89638e-07 ,6.14115e-07 ,6.39320e-07 , &
      6.65266e-07 ,6.91965e-07 ,7.19427e-07 ,7.47666e-07 ,7.76691e-07 , &
      8.06516e-07 ,8.37151e-07 ,8.68607e-07 ,9.00896e-07 ,9.34029e-07 , &
      9.68018e-07 ,1.00287e-06 ,1.03860e-06 ,1.07522e-06 ,1.11274e-06 , &
      1.15117e-06 ,1.19052e-06 ,1.23079e-06 ,1.27201e-06 ,1.31418e-06 , &
      1.35731e-06 ,1.40141e-06 ,1.44650e-06 ,1.49257e-06 ,1.53965e-06 , &
      1.58773e-06 ,1.63684e-06 ,1.68697e-06 ,1.73815e-06 ,1.79037e-06 /)
      totplnk(51:100, 6) = (/ &
      1.84365e-06 ,1.89799e-06 ,1.95341e-06 ,2.00991e-06 ,2.06750e-06 , &
      2.12619e-06 ,2.18599e-06 ,2.24691e-06 ,2.30895e-06 ,2.37212e-06 , &
      2.43643e-06 ,2.50189e-06 ,2.56851e-06 ,2.63628e-06 ,2.70523e-06 , &
      2.77536e-06 ,2.84666e-06 ,2.91916e-06 ,2.99286e-06 ,3.06776e-06 , &
      3.14387e-06 ,3.22120e-06 ,3.29975e-06 ,3.37953e-06 ,3.46054e-06 , &
      3.54280e-06 ,3.62630e-06 ,3.71105e-06 ,3.79707e-06 ,3.88434e-06 , &
      3.97288e-06 ,4.06270e-06 ,4.15380e-06 ,4.24617e-06 ,4.33984e-06 , &
      4.43479e-06 ,4.53104e-06 ,4.62860e-06 ,4.72746e-06 ,4.82763e-06 , &
      4.92911e-06 ,5.03191e-06 ,5.13603e-06 ,5.24147e-06 ,5.34824e-06 , &
      5.45634e-06 ,5.56578e-06 ,5.67656e-06 ,5.78867e-06 ,5.90213e-06 /)
      totplnk(101:150, 6) = (/ &
      6.01694e-06 ,6.13309e-06 ,6.25060e-06 ,6.36947e-06 ,6.48968e-06 , &
      6.61126e-06 ,6.73420e-06 ,6.85850e-06 ,6.98417e-06 ,7.11120e-06 , &
      7.23961e-06 ,7.36938e-06 ,7.50053e-06 ,7.63305e-06 ,7.76694e-06 , &
      7.90221e-06 ,8.03887e-06 ,8.17690e-06 ,8.31632e-06 ,8.45710e-06 , &
      8.59928e-06 ,8.74282e-06 ,8.88776e-06 ,9.03409e-06 ,9.18179e-06 , &
      9.33088e-06 ,9.48136e-06 ,9.63323e-06 ,9.78648e-06 ,9.94111e-06 , &
      1.00971e-05 ,1.02545e-05 ,1.04133e-05 ,1.05735e-05 ,1.07351e-05 , &
      1.08980e-05 ,1.10624e-05 ,1.12281e-05 ,1.13952e-05 ,1.15637e-05 , &
      1.17335e-05 ,1.19048e-05 ,1.20774e-05 ,1.22514e-05 ,1.24268e-05 , &
      1.26036e-05 ,1.27817e-05 ,1.29612e-05 ,1.31421e-05 ,1.33244e-05 /)
      totplnk(151:181, 6) = (/ &
      1.35080e-05 ,1.36930e-05 ,1.38794e-05 ,1.40672e-05 ,1.42563e-05 , &
      1.44468e-05 ,1.46386e-05 ,1.48318e-05 ,1.50264e-05 ,1.52223e-05 , &
      1.54196e-05 ,1.56182e-05 ,1.58182e-05 ,1.60196e-05 ,1.62223e-05 , &
      1.64263e-05 ,1.66317e-05 ,1.68384e-05 ,1.70465e-05 ,1.72559e-05 , &
      1.74666e-05 ,1.76787e-05 ,1.78921e-05 ,1.81069e-05 ,1.83230e-05 , &
      1.85404e-05 ,1.87591e-05 ,1.89791e-05 ,1.92005e-05 ,1.94232e-05 , &
      1.96471e-05 /)
      totplnk(1:50, 7) = (/ &
      1.25349e-07 ,1.32735e-07 ,1.40458e-07 ,1.48527e-07 ,1.56954e-07 , &
      1.65748e-07 ,1.74920e-07 ,1.84481e-07 ,1.94443e-07 ,2.04814e-07 , &
      2.15608e-07 ,2.26835e-07 ,2.38507e-07 ,2.50634e-07 ,2.63229e-07 , &
      2.76301e-07 ,2.89864e-07 ,3.03930e-07 ,3.18508e-07 ,3.33612e-07 , &
      3.49253e-07 ,3.65443e-07 ,3.82195e-07 ,3.99519e-07 ,4.17428e-07 , &
      4.35934e-07 ,4.55050e-07 ,4.74785e-07 ,4.95155e-07 ,5.16170e-07 , &
      5.37844e-07 ,5.60186e-07 ,5.83211e-07 ,6.06929e-07 ,6.31355e-07 , &
      6.56498e-07 ,6.82373e-07 ,7.08990e-07 ,7.36362e-07 ,7.64501e-07 , &
      7.93420e-07 ,8.23130e-07 ,8.53643e-07 ,8.84971e-07 ,9.17128e-07 , &
      9.50123e-07 ,9.83969e-07 ,1.01868e-06 ,1.05426e-06 ,1.09073e-06 /)
      totplnk(51:100, 7) = (/ &
      1.12810e-06 ,1.16638e-06 ,1.20558e-06 ,1.24572e-06 ,1.28680e-06 , &
      1.32883e-06 ,1.37183e-06 ,1.41581e-06 ,1.46078e-06 ,1.50675e-06 , &
      1.55374e-06 ,1.60174e-06 ,1.65078e-06 ,1.70087e-06 ,1.75200e-06 , &
      1.80421e-06 ,1.85749e-06 ,1.91186e-06 ,1.96732e-06 ,2.02389e-06 , &
      2.08159e-06 ,2.14040e-06 ,2.20035e-06 ,2.26146e-06 ,2.32372e-06 , &
      2.38714e-06 ,2.45174e-06 ,2.51753e-06 ,2.58451e-06 ,2.65270e-06 , &
      2.72210e-06 ,2.79272e-06 ,2.86457e-06 ,2.93767e-06 ,3.01201e-06 , &
      3.08761e-06 ,3.16448e-06 ,3.24261e-06 ,3.32204e-06 ,3.40275e-06 , &
      3.48476e-06 ,3.56808e-06 ,3.65271e-06 ,3.73866e-06 ,3.82595e-06 , &
      3.91456e-06 ,4.00453e-06 ,4.09584e-06 ,4.18851e-06 ,4.28254e-06 /)
      totplnk(101:150, 7) = (/ &
      4.37796e-06 ,4.47475e-06 ,4.57293e-06 ,4.67249e-06 ,4.77346e-06 , &
      4.87583e-06 ,4.97961e-06 ,5.08481e-06 ,5.19143e-06 ,5.29948e-06 , &
      5.40896e-06 ,5.51989e-06 ,5.63226e-06 ,5.74608e-06 ,5.86136e-06 , &
      5.97810e-06 ,6.09631e-06 ,6.21597e-06 ,6.33713e-06 ,6.45976e-06 , &
      6.58388e-06 ,6.70950e-06 ,6.83661e-06 ,6.96521e-06 ,7.09531e-06 , &
      7.22692e-06 ,7.36005e-06 ,7.49468e-06 ,7.63084e-06 ,7.76851e-06 , &
      7.90773e-06 ,8.04846e-06 ,8.19072e-06 ,8.33452e-06 ,8.47985e-06 , &
      8.62674e-06 ,8.77517e-06 ,8.92514e-06 ,9.07666e-06 ,9.22975e-06 , &
      9.38437e-06 ,9.54057e-06 ,9.69832e-06 ,9.85762e-06 ,1.00185e-05 , &
      1.01810e-05 ,1.03450e-05 ,1.05106e-05 ,1.06777e-05 ,1.08465e-05 /)
      totplnk(151:181, 7) = (/ &
      1.10168e-05 ,1.11887e-05 ,1.13621e-05 ,1.15372e-05 ,1.17138e-05 , &
      1.18920e-05 ,1.20718e-05 ,1.22532e-05 ,1.24362e-05 ,1.26207e-05 , &
      1.28069e-05 ,1.29946e-05 ,1.31839e-05 ,1.33749e-05 ,1.35674e-05 , &
      1.37615e-05 ,1.39572e-05 ,1.41544e-05 ,1.43533e-05 ,1.45538e-05 , &
      1.47558e-05 ,1.49595e-05 ,1.51647e-05 ,1.53716e-05 ,1.55800e-05 , &
      1.57900e-05 ,1.60017e-05 ,1.62149e-05 ,1.64296e-05 ,1.66460e-05 , &
      1.68640e-05 /)
      totplnk(1:50, 8) = (/ &
      6.74445e-08 ,7.18176e-08 ,7.64153e-08 ,8.12456e-08 ,8.63170e-08 , &
      9.16378e-08 ,9.72168e-08 ,1.03063e-07 ,1.09184e-07 ,1.15591e-07 , &
      1.22292e-07 ,1.29296e-07 ,1.36613e-07 ,1.44253e-07 ,1.52226e-07 , &
      1.60540e-07 ,1.69207e-07 ,1.78236e-07 ,1.87637e-07 ,1.97421e-07 , &
      2.07599e-07 ,2.18181e-07 ,2.29177e-07 ,2.40598e-07 ,2.52456e-07 , &
      2.64761e-07 ,2.77523e-07 ,2.90755e-07 ,3.04468e-07 ,3.18673e-07 , &
      3.33381e-07 ,3.48603e-07 ,3.64352e-07 ,3.80638e-07 ,3.97474e-07 , &
      4.14871e-07 ,4.32841e-07 ,4.51395e-07 ,4.70547e-07 ,4.90306e-07 , &
      5.10687e-07 ,5.31699e-07 ,5.53357e-07 ,5.75670e-07 ,5.98652e-07 , &
      6.22315e-07 ,6.46672e-07 ,6.71731e-07 ,6.97511e-07 ,7.24018e-07 /)
      totplnk(51:100, 8) = (/ &
      7.51266e-07 ,7.79269e-07 ,8.08038e-07 ,8.37584e-07 ,8.67922e-07 , &
      8.99061e-07 ,9.31016e-07 ,9.63797e-07 ,9.97417e-07 ,1.03189e-06 , &
      1.06722e-06 ,1.10343e-06 ,1.14053e-06 ,1.17853e-06 ,1.21743e-06 , &
      1.25726e-06 ,1.29803e-06 ,1.33974e-06 ,1.38241e-06 ,1.42606e-06 , &
      1.47068e-06 ,1.51630e-06 ,1.56293e-06 ,1.61056e-06 ,1.65924e-06 , &
      1.70894e-06 ,1.75971e-06 ,1.81153e-06 ,1.86443e-06 ,1.91841e-06 , &
      1.97350e-06 ,2.02968e-06 ,2.08699e-06 ,2.14543e-06 ,2.20500e-06 , &
      2.26573e-06 ,2.32762e-06 ,2.39068e-06 ,2.45492e-06 ,2.52036e-06 , &
      2.58700e-06 ,2.65485e-06 ,2.72393e-06 ,2.79424e-06 ,2.86580e-06 , &
      2.93861e-06 ,3.01269e-06 ,3.08803e-06 ,3.16467e-06 ,3.24259e-06 /)
      totplnk(101:150, 8) = (/ &
      3.32181e-06 ,3.40235e-06 ,3.48420e-06 ,3.56739e-06 ,3.65192e-06 , &
      3.73779e-06 ,3.82502e-06 ,3.91362e-06 ,4.00359e-06 ,4.09494e-06 , &
      4.18768e-06 ,4.28182e-06 ,4.37737e-06 ,4.47434e-06 ,4.57273e-06 , &
      4.67254e-06 ,4.77380e-06 ,4.87651e-06 ,4.98067e-06 ,5.08630e-06 , &
      5.19339e-06 ,5.30196e-06 ,5.41201e-06 ,5.52356e-06 ,5.63660e-06 , &
      5.75116e-06 ,5.86722e-06 ,5.98479e-06 ,6.10390e-06 ,6.22453e-06 , &
      6.34669e-06 ,6.47042e-06 ,6.59569e-06 ,6.72252e-06 ,6.85090e-06 , &
      6.98085e-06 ,7.11238e-06 ,7.24549e-06 ,7.38019e-06 ,7.51646e-06 , &
      7.65434e-06 ,7.79382e-06 ,7.93490e-06 ,8.07760e-06 ,8.22192e-06 , &
      8.36784e-06 ,8.51540e-06 ,8.66459e-06 ,8.81542e-06 ,8.96786e-06 /)
      totplnk(151:181, 8) = (/ &
      9.12197e-06 ,9.27772e-06 ,9.43513e-06 ,9.59419e-06 ,9.75490e-06 , &
      9.91728e-06 ,1.00813e-05 ,1.02471e-05 ,1.04144e-05 ,1.05835e-05 , &
      1.07543e-05 ,1.09267e-05 ,1.11008e-05 ,1.12766e-05 ,1.14541e-05 , &
      1.16333e-05 ,1.18142e-05 ,1.19969e-05 ,1.21812e-05 ,1.23672e-05 , &
      1.25549e-05 ,1.27443e-05 ,1.29355e-05 ,1.31284e-05 ,1.33229e-05 , &
      1.35193e-05 ,1.37173e-05 ,1.39170e-05 ,1.41185e-05 ,1.43217e-05 , &
      1.45267e-05 /)
      totplnk(1:50, 9) = (/ &
      2.61522e-08 ,2.80613e-08 ,3.00838e-08 ,3.22250e-08 ,3.44899e-08 , &
      3.68841e-08 ,3.94129e-08 ,4.20820e-08 ,4.48973e-08 ,4.78646e-08 , &
      5.09901e-08 ,5.42799e-08 ,5.77405e-08 ,6.13784e-08 ,6.52001e-08 , &
      6.92126e-08 ,7.34227e-08 ,7.78375e-08 ,8.24643e-08 ,8.73103e-08 , &
      9.23832e-08 ,9.76905e-08 ,1.03240e-07 ,1.09039e-07 ,1.15097e-07 , &
      1.21421e-07 ,1.28020e-07 ,1.34902e-07 ,1.42075e-07 ,1.49548e-07 , &
      1.57331e-07 ,1.65432e-07 ,1.73860e-07 ,1.82624e-07 ,1.91734e-07 , &
      2.01198e-07 ,2.11028e-07 ,2.21231e-07 ,2.31818e-07 ,2.42799e-07 , &
      2.54184e-07 ,2.65983e-07 ,2.78205e-07 ,2.90862e-07 ,3.03963e-07 , &
      3.17519e-07 ,3.31541e-07 ,3.46039e-07 ,3.61024e-07 ,3.76507e-07 /)
      totplnk(51:100, 9) = (/ &
      3.92498e-07 ,4.09008e-07 ,4.26050e-07 ,4.43633e-07 ,4.61769e-07 , &
      4.80469e-07 ,4.99744e-07 ,5.19606e-07 ,5.40067e-07 ,5.61136e-07 , &
      5.82828e-07 ,6.05152e-07 ,6.28120e-07 ,6.51745e-07 ,6.76038e-07 , &
      7.01010e-07 ,7.26674e-07 ,7.53041e-07 ,7.80124e-07 ,8.07933e-07 , &
      8.36482e-07 ,8.65781e-07 ,8.95845e-07 ,9.26683e-07 ,9.58308e-07 , &
      9.90732e-07 ,1.02397e-06 ,1.05803e-06 ,1.09292e-06 ,1.12866e-06 , &
      1.16526e-06 ,1.20274e-06 ,1.24109e-06 ,1.28034e-06 ,1.32050e-06 , &
      1.36158e-06 ,1.40359e-06 ,1.44655e-06 ,1.49046e-06 ,1.53534e-06 , &
      1.58120e-06 ,1.62805e-06 ,1.67591e-06 ,1.72478e-06 ,1.77468e-06 , &
      1.82561e-06 ,1.87760e-06 ,1.93066e-06 ,1.98479e-06 ,2.04000e-06 /)
      totplnk(101:150, 9) = (/ &
      2.09631e-06 ,2.15373e-06 ,2.21228e-06 ,2.27196e-06 ,2.33278e-06 , &
      2.39475e-06 ,2.45790e-06 ,2.52222e-06 ,2.58773e-06 ,2.65445e-06 , &
      2.72238e-06 ,2.79152e-06 ,2.86191e-06 ,2.93354e-06 ,3.00643e-06 , &
      3.08058e-06 ,3.15601e-06 ,3.23273e-06 ,3.31075e-06 ,3.39009e-06 , &
      3.47074e-06 ,3.55272e-06 ,3.63605e-06 ,3.72072e-06 ,3.80676e-06 , &
      3.89417e-06 ,3.98297e-06 ,4.07315e-06 ,4.16474e-06 ,4.25774e-06 , &
      4.35217e-06 ,4.44802e-06 ,4.54532e-06 ,4.64406e-06 ,4.74428e-06 , &
      4.84595e-06 ,4.94911e-06 ,5.05376e-06 ,5.15990e-06 ,5.26755e-06 , &
      5.37671e-06 ,5.48741e-06 ,5.59963e-06 ,5.71340e-06 ,5.82871e-06 , &
      5.94559e-06 ,6.06403e-06 ,6.18404e-06 ,6.30565e-06 ,6.42885e-06 /)
      totplnk(151:181, 9) = (/ &
      6.55364e-06 ,6.68004e-06 ,6.80806e-06 ,6.93771e-06 ,7.06898e-06 , &
      7.20190e-06 ,7.33646e-06 ,7.47267e-06 ,7.61056e-06 ,7.75010e-06 , &
      7.89133e-06 ,8.03423e-06 ,8.17884e-06 ,8.32514e-06 ,8.47314e-06 , &
      8.62284e-06 ,8.77427e-06 ,8.92743e-06 ,9.08231e-06 ,9.23893e-06 , &
      9.39729e-06 ,9.55741e-06 ,9.71927e-06 ,9.88291e-06 ,1.00483e-05 , &
      1.02155e-05 ,1.03844e-05 ,1.05552e-05 ,1.07277e-05 ,1.09020e-05 , &
      1.10781e-05 /)
      totplnk(1:50,10) = (/ &
      8.89300e-09 ,9.63263e-09 ,1.04235e-08 ,1.12685e-08 ,1.21703e-08 , &
      1.31321e-08 ,1.41570e-08 ,1.52482e-08 ,1.64090e-08 ,1.76428e-08 , &
      1.89533e-08 ,2.03441e-08 ,2.18190e-08 ,2.33820e-08 ,2.50370e-08 , &
      2.67884e-08 ,2.86402e-08 ,3.05969e-08 ,3.26632e-08 ,3.48436e-08 , &
      3.71429e-08 ,3.95660e-08 ,4.21179e-08 ,4.48040e-08 ,4.76294e-08 , &
      5.05996e-08 ,5.37201e-08 ,5.69966e-08 ,6.04349e-08 ,6.40411e-08 , &
      6.78211e-08 ,7.17812e-08 ,7.59276e-08 ,8.02670e-08 ,8.48059e-08 , &
      8.95508e-08 ,9.45090e-08 ,9.96873e-08 ,1.05093e-07 ,1.10733e-07 , &
      1.16614e-07 ,1.22745e-07 ,1.29133e-07 ,1.35786e-07 ,1.42711e-07 , &
      1.49916e-07 ,1.57410e-07 ,1.65202e-07 ,1.73298e-07 ,1.81709e-07 /)
      totplnk(51:100,10) = (/ &
      1.90441e-07 ,1.99505e-07 ,2.08908e-07 ,2.18660e-07 ,2.28770e-07 , &
      2.39247e-07 ,2.50101e-07 ,2.61340e-07 ,2.72974e-07 ,2.85013e-07 , &
      2.97467e-07 ,3.10345e-07 ,3.23657e-07 ,3.37413e-07 ,3.51623e-07 , &
      3.66298e-07 ,3.81448e-07 ,3.97082e-07 ,4.13212e-07 ,4.29848e-07 , &
      4.47000e-07 ,4.64680e-07 ,4.82898e-07 ,5.01664e-07 ,5.20991e-07 , &
      5.40888e-07 ,5.61369e-07 ,5.82440e-07 ,6.04118e-07 ,6.26410e-07 , &
      6.49329e-07 ,6.72887e-07 ,6.97095e-07 ,7.21964e-07 ,7.47506e-07 , &
      7.73732e-07 ,8.00655e-07 ,8.28287e-07 ,8.56635e-07 ,8.85717e-07 , &
      9.15542e-07 ,9.46122e-07 ,9.77469e-07 ,1.00960e-06 ,1.04251e-06 , &
      1.07623e-06 ,1.11077e-06 ,1.14613e-06 ,1.18233e-06 ,1.21939e-06 /)
      totplnk(101:150,10) = (/ &
      1.25730e-06 ,1.29610e-06 ,1.33578e-06 ,1.37636e-06 ,1.41785e-06 , &
      1.46027e-06 ,1.50362e-06 ,1.54792e-06 ,1.59319e-06 ,1.63942e-06 , &
      1.68665e-06 ,1.73487e-06 ,1.78410e-06 ,1.83435e-06 ,1.88564e-06 , &
      1.93797e-06 ,1.99136e-06 ,2.04582e-06 ,2.10137e-06 ,2.15801e-06 , &
      2.21576e-06 ,2.27463e-06 ,2.33462e-06 ,2.39577e-06 ,2.45806e-06 , &
      2.52153e-06 ,2.58617e-06 ,2.65201e-06 ,2.71905e-06 ,2.78730e-06 , &
      2.85678e-06 ,2.92749e-06 ,2.99946e-06 ,3.07269e-06 ,3.14720e-06 , &
      3.22299e-06 ,3.30007e-06 ,3.37847e-06 ,3.45818e-06 ,3.53923e-06 , &
      3.62161e-06 ,3.70535e-06 ,3.79046e-06 ,3.87695e-06 ,3.96481e-06 , &
      4.05409e-06 ,4.14477e-06 ,4.23687e-06 ,4.33040e-06 ,4.42538e-06 /)
      totplnk(151:181,10) = (/ &
      4.52180e-06 ,4.61969e-06 ,4.71905e-06 ,4.81991e-06 ,4.92226e-06 , &
      5.02611e-06 ,5.13148e-06 ,5.23839e-06 ,5.34681e-06 ,5.45681e-06 , &
      5.56835e-06 ,5.68146e-06 ,5.79614e-06 ,5.91242e-06 ,6.03030e-06 , &
      6.14978e-06 ,6.27088e-06 ,6.39360e-06 ,6.51798e-06 ,6.64398e-06 , &
      6.77165e-06 ,6.90099e-06 ,7.03198e-06 ,7.16468e-06 ,7.29906e-06 , &
      7.43514e-06 ,7.57294e-06 ,7.71244e-06 ,7.85369e-06 ,7.99666e-06 , &
      8.14138e-06 /)
      totplnk(1:50,11) = (/ &
      2.53767e-09 ,2.77242e-09 ,3.02564e-09 ,3.29851e-09 ,3.59228e-09 , &
      3.90825e-09 ,4.24777e-09 ,4.61227e-09 ,5.00322e-09 ,5.42219e-09 , &
      5.87080e-09 ,6.35072e-09 ,6.86370e-09 ,7.41159e-09 ,7.99628e-09 , &
      8.61974e-09 ,9.28404e-09 ,9.99130e-09 ,1.07437e-08 ,1.15436e-08 , &
      1.23933e-08 ,1.32953e-08 ,1.42522e-08 ,1.52665e-08 ,1.63410e-08 , &
      1.74786e-08 ,1.86820e-08 ,1.99542e-08 ,2.12985e-08 ,2.27179e-08 , &
      2.42158e-08 ,2.57954e-08 ,2.74604e-08 ,2.92141e-08 ,3.10604e-08 , &
      3.30029e-08 ,3.50457e-08 ,3.71925e-08 ,3.94476e-08 ,4.18149e-08 , &
      4.42991e-08 ,4.69043e-08 ,4.96352e-08 ,5.24961e-08 ,5.54921e-08 , &
      5.86277e-08 ,6.19081e-08 ,6.53381e-08 ,6.89231e-08 ,7.26681e-08 /)
      totplnk(51:100,11) = (/ &
      7.65788e-08 ,8.06604e-08 ,8.49187e-08 ,8.93591e-08 ,9.39879e-08 , &
      9.88106e-08 ,1.03834e-07 ,1.09063e-07 ,1.14504e-07 ,1.20165e-07 , &
      1.26051e-07 ,1.32169e-07 ,1.38525e-07 ,1.45128e-07 ,1.51982e-07 , &
      1.59096e-07 ,1.66477e-07 ,1.74132e-07 ,1.82068e-07 ,1.90292e-07 , &
      1.98813e-07 ,2.07638e-07 ,2.16775e-07 ,2.26231e-07 ,2.36015e-07 , &
      2.46135e-07 ,2.56599e-07 ,2.67415e-07 ,2.78592e-07 ,2.90137e-07 , &
      3.02061e-07 ,3.14371e-07 ,3.27077e-07 ,3.40186e-07 ,3.53710e-07 , &
      3.67655e-07 ,3.82031e-07 ,3.96848e-07 ,4.12116e-07 ,4.27842e-07 , &
      4.44039e-07 ,4.60713e-07 ,4.77876e-07 ,4.95537e-07 ,5.13706e-07 , &
      5.32392e-07 ,5.51608e-07 ,5.71360e-07 ,5.91662e-07 ,6.12521e-07 /)
      totplnk(101:150,11) = (/ &
      6.33950e-07 ,6.55958e-07 ,6.78556e-07 ,7.01753e-07 ,7.25562e-07 , &
      7.49992e-07 ,7.75055e-07 ,8.00760e-07 ,8.27120e-07 ,8.54145e-07 , &
      8.81845e-07 ,9.10233e-07 ,9.39318e-07 ,9.69113e-07 ,9.99627e-07 , &
      1.03087e-06 ,1.06286e-06 ,1.09561e-06 ,1.12912e-06 ,1.16340e-06 , &
      1.19848e-06 ,1.23435e-06 ,1.27104e-06 ,1.30855e-06 ,1.34690e-06 , &
      1.38609e-06 ,1.42614e-06 ,1.46706e-06 ,1.50886e-06 ,1.55155e-06 , &
      1.59515e-06 ,1.63967e-06 ,1.68512e-06 ,1.73150e-06 ,1.77884e-06 , &
      1.82715e-06 ,1.87643e-06 ,1.92670e-06 ,1.97797e-06 ,2.03026e-06 , &
      2.08356e-06 ,2.13791e-06 ,2.19330e-06 ,2.24975e-06 ,2.30728e-06 , &
      2.36589e-06 ,2.42560e-06 ,2.48641e-06 ,2.54835e-06 ,2.61142e-06 /)
      totplnk(151:181,11) = (/ &
      2.67563e-06 ,2.74100e-06 ,2.80754e-06 ,2.87526e-06 ,2.94417e-06 , &
      3.01429e-06 ,3.08562e-06 ,3.15819e-06 ,3.23199e-06 ,3.30704e-06 , &
      3.38336e-06 ,3.46096e-06 ,3.53984e-06 ,3.62002e-06 ,3.70151e-06 , &
      3.78433e-06 ,3.86848e-06 ,3.95399e-06 ,4.04084e-06 ,4.12907e-06 , &
      4.21868e-06 ,4.30968e-06 ,4.40209e-06 ,4.49592e-06 ,4.59117e-06 , &
      4.68786e-06 ,4.78600e-06 ,4.88561e-06 ,4.98669e-06 ,5.08926e-06 , &
      5.19332e-06 /)
      totplnk(1:50,12) = (/ &
      2.73921e-10 ,3.04500e-10 ,3.38056e-10 ,3.74835e-10 ,4.15099e-10 , &
      4.59126e-10 ,5.07214e-10 ,5.59679e-10 ,6.16857e-10 ,6.79103e-10 , &
      7.46796e-10 ,8.20335e-10 ,9.00144e-10 ,9.86671e-10 ,1.08039e-09 , &
      1.18180e-09 ,1.29142e-09 ,1.40982e-09 ,1.53757e-09 ,1.67529e-09 , &
      1.82363e-09 ,1.98327e-09 ,2.15492e-09 ,2.33932e-09 ,2.53726e-09 , &
      2.74957e-09 ,2.97710e-09 ,3.22075e-09 ,3.48145e-09 ,3.76020e-09 , &
      4.05801e-09 ,4.37595e-09 ,4.71513e-09 ,5.07672e-09 ,5.46193e-09 , &
      5.87201e-09 ,6.30827e-09 ,6.77205e-09 ,7.26480e-09 ,7.78794e-09 , &
      8.34304e-09 ,8.93163e-09 ,9.55537e-09 ,1.02159e-08 ,1.09151e-08 , &
      1.16547e-08 ,1.24365e-08 ,1.32625e-08 ,1.41348e-08 ,1.50554e-08 /)
      totplnk(51:100,12) = (/ &
      1.60264e-08 ,1.70500e-08 ,1.81285e-08 ,1.92642e-08 ,2.04596e-08 , &
      2.17171e-08 ,2.30394e-08 ,2.44289e-08 ,2.58885e-08 ,2.74209e-08 , &
      2.90290e-08 ,3.07157e-08 ,3.24841e-08 ,3.43371e-08 ,3.62782e-08 , &
      3.83103e-08 ,4.04371e-08 ,4.26617e-08 ,4.49878e-08 ,4.74190e-08 , &
      4.99589e-08 ,5.26113e-08 ,5.53801e-08 ,5.82692e-08 ,6.12826e-08 , &
      6.44245e-08 ,6.76991e-08 ,7.11105e-08 ,7.46634e-08 ,7.83621e-08 , &
      8.22112e-08 ,8.62154e-08 ,9.03795e-08 ,9.47081e-08 ,9.92066e-08 , &
      1.03879e-07 ,1.08732e-07 ,1.13770e-07 ,1.18998e-07 ,1.24422e-07 , &
      1.30048e-07 ,1.35880e-07 ,1.41924e-07 ,1.48187e-07 ,1.54675e-07 , &
      1.61392e-07 ,1.68346e-07 ,1.75543e-07 ,1.82988e-07 ,1.90688e-07 /)
      totplnk(101:150,12) = (/ &
      1.98650e-07 ,2.06880e-07 ,2.15385e-07 ,2.24172e-07 ,2.33247e-07 , &
      2.42617e-07 ,2.52289e-07 ,2.62272e-07 ,2.72571e-07 ,2.83193e-07 , &
      2.94147e-07 ,3.05440e-07 ,3.17080e-07 ,3.29074e-07 ,3.41430e-07 , &
      3.54155e-07 ,3.67259e-07 ,3.80747e-07 ,3.94631e-07 ,4.08916e-07 , &
      4.23611e-07 ,4.38725e-07 ,4.54267e-07 ,4.70245e-07 ,4.86666e-07 , &
      5.03541e-07 ,5.20879e-07 ,5.38687e-07 ,5.56975e-07 ,5.75751e-07 , &
      5.95026e-07 ,6.14808e-07 ,6.35107e-07 ,6.55932e-07 ,6.77293e-07 , &
      6.99197e-07 ,7.21656e-07 ,7.44681e-07 ,7.68278e-07 ,7.92460e-07 , &
      8.17235e-07 ,8.42614e-07 ,8.68606e-07 ,8.95223e-07 ,9.22473e-07 , &
      9.50366e-07 ,9.78915e-07 ,1.00813e-06 ,1.03802e-06 ,1.06859e-06 /)
      totplnk(151:181,12) = (/ &
      1.09986e-06 ,1.13184e-06 ,1.16453e-06 ,1.19796e-06 ,1.23212e-06 , &
      1.26703e-06 ,1.30270e-06 ,1.33915e-06 ,1.37637e-06 ,1.41440e-06 , &
      1.45322e-06 ,1.49286e-06 ,1.53333e-06 ,1.57464e-06 ,1.61679e-06 , &
      1.65981e-06 ,1.70370e-06 ,1.74847e-06 ,1.79414e-06 ,1.84071e-06 , &
      1.88821e-06 ,1.93663e-06 ,1.98599e-06 ,2.03631e-06 ,2.08759e-06 , &
      2.13985e-06 ,2.19310e-06 ,2.24734e-06 ,2.30260e-06 ,2.35888e-06 , &
      2.41619e-06 /)
      totplnk(1:50,13) = (/ &
      4.53634e-11 ,5.11435e-11 ,5.75754e-11 ,6.47222e-11 ,7.26531e-11 , &
      8.14420e-11 ,9.11690e-11 ,1.01921e-10 ,1.13790e-10 ,1.26877e-10 , &
      1.41288e-10 ,1.57140e-10 ,1.74555e-10 ,1.93665e-10 ,2.14613e-10 , &
      2.37548e-10 ,2.62633e-10 ,2.90039e-10 ,3.19948e-10 ,3.52558e-10 , &
      3.88073e-10 ,4.26716e-10 ,4.68719e-10 ,5.14331e-10 ,5.63815e-10 , &
      6.17448e-10 ,6.75526e-10 ,7.38358e-10 ,8.06277e-10 ,8.79625e-10 , &
      9.58770e-10 ,1.04410e-09 ,1.13602e-09 ,1.23495e-09 ,1.34135e-09 , &
      1.45568e-09 ,1.57845e-09 ,1.71017e-09 ,1.85139e-09 ,2.00268e-09 , &
      2.16464e-09 ,2.33789e-09 ,2.52309e-09 ,2.72093e-09 ,2.93212e-09 , &
      3.15740e-09 ,3.39757e-09 ,3.65341e-09 ,3.92579e-09 ,4.21559e-09 /)
      totplnk(51:100,13) = (/ &
      4.52372e-09 ,4.85115e-09 ,5.19886e-09 ,5.56788e-09 ,5.95928e-09 , &
      6.37419e-09 ,6.81375e-09 ,7.27917e-09 ,7.77168e-09 ,8.29256e-09 , &
      8.84317e-09 ,9.42487e-09 ,1.00391e-08 ,1.06873e-08 ,1.13710e-08 , &
      1.20919e-08 ,1.28515e-08 ,1.36514e-08 ,1.44935e-08 ,1.53796e-08 , &
      1.63114e-08 ,1.72909e-08 ,1.83201e-08 ,1.94008e-08 ,2.05354e-08 , &
      2.17258e-08 ,2.29742e-08 ,2.42830e-08 ,2.56545e-08 ,2.70910e-08 , &
      2.85950e-08 ,3.01689e-08 ,3.18155e-08 ,3.35373e-08 ,3.53372e-08 , &
      3.72177e-08 ,3.91818e-08 ,4.12325e-08 ,4.33727e-08 ,4.56056e-08 , &
      4.79342e-08 ,5.03617e-08 ,5.28915e-08 ,5.55270e-08 ,5.82715e-08 , &
      6.11286e-08 ,6.41019e-08 ,6.71951e-08 ,7.04119e-08 ,7.37560e-08 /)
      totplnk(101:150,13) = (/ &
      7.72315e-08 ,8.08424e-08 ,8.45927e-08 ,8.84866e-08 ,9.25281e-08 , &
      9.67218e-08 ,1.01072e-07 ,1.05583e-07 ,1.10260e-07 ,1.15107e-07 , &
      1.20128e-07 ,1.25330e-07 ,1.30716e-07 ,1.36291e-07 ,1.42061e-07 , &
      1.48031e-07 ,1.54206e-07 ,1.60592e-07 ,1.67192e-07 ,1.74015e-07 , &
      1.81064e-07 ,1.88345e-07 ,1.95865e-07 ,2.03628e-07 ,2.11643e-07 , &
      2.19912e-07 ,2.28443e-07 ,2.37244e-07 ,2.46318e-07 ,2.55673e-07 , &
      2.65316e-07 ,2.75252e-07 ,2.85489e-07 ,2.96033e-07 ,3.06891e-07 , &
      3.18070e-07 ,3.29576e-07 ,3.41417e-07 ,3.53600e-07 ,3.66133e-07 , &
      3.79021e-07 ,3.92274e-07 ,4.05897e-07 ,4.19899e-07 ,4.34288e-07 , &
      4.49071e-07 ,4.64255e-07 ,4.79850e-07 ,4.95863e-07 ,5.12300e-07 /)
      totplnk(151:181,13) = (/ &
      5.29172e-07 ,5.46486e-07 ,5.64250e-07 ,5.82473e-07 ,6.01164e-07 , &
      6.20329e-07 ,6.39979e-07 ,6.60122e-07 ,6.80767e-07 ,7.01922e-07 , &
      7.23596e-07 ,7.45800e-07 ,7.68539e-07 ,7.91826e-07 ,8.15669e-07 , &
      8.40076e-07 ,8.65058e-07 ,8.90623e-07 ,9.16783e-07 ,9.43544e-07 , &
      9.70917e-07 ,9.98912e-07 ,1.02754e-06 ,1.05681e-06 ,1.08673e-06 , &
      1.11731e-06 ,1.14856e-06 ,1.18050e-06 ,1.21312e-06 ,1.24645e-06 , &
      1.28049e-06 /)
      totplnk(1:50,14) = (/ &
      1.40113e-11 ,1.59358e-11 ,1.80960e-11 ,2.05171e-11 ,2.32266e-11 , &
      2.62546e-11 ,2.96335e-11 ,3.33990e-11 ,3.75896e-11 ,4.22469e-11 , &
      4.74164e-11 ,5.31466e-11 ,5.94905e-11 ,6.65054e-11 ,7.42522e-11 , &
      8.27975e-11 ,9.22122e-11 ,1.02573e-10 ,1.13961e-10 ,1.26466e-10 , &
      1.40181e-10 ,1.55206e-10 ,1.71651e-10 ,1.89630e-10 ,2.09265e-10 , &
      2.30689e-10 ,2.54040e-10 ,2.79467e-10 ,3.07128e-10 ,3.37190e-10 , &
      3.69833e-10 ,4.05243e-10 ,4.43623e-10 ,4.85183e-10 ,5.30149e-10 , &
      5.78755e-10 ,6.31255e-10 ,6.87910e-10 ,7.49002e-10 ,8.14824e-10 , &
      8.85687e-10 ,9.61914e-10 ,1.04385e-09 ,1.13186e-09 ,1.22631e-09 , &
      1.32761e-09 ,1.43617e-09 ,1.55243e-09 ,1.67686e-09 ,1.80992e-09 /)
      totplnk(51:100,14) = (/ &
      1.95212e-09 ,2.10399e-09 ,2.26607e-09 ,2.43895e-09 ,2.62321e-09 , &
      2.81949e-09 ,3.02844e-09 ,3.25073e-09 ,3.48707e-09 ,3.73820e-09 , &
      4.00490e-09 ,4.28794e-09 ,4.58819e-09 ,4.90647e-09 ,5.24371e-09 , &
      5.60081e-09 ,5.97875e-09 ,6.37854e-09 ,6.80120e-09 ,7.24782e-09 , &
      7.71950e-09 ,8.21740e-09 ,8.74271e-09 ,9.29666e-09 ,9.88054e-09 , &
      1.04956e-08 ,1.11434e-08 ,1.18251e-08 ,1.25422e-08 ,1.32964e-08 , &
      1.40890e-08 ,1.49217e-08 ,1.57961e-08 ,1.67140e-08 ,1.76771e-08 , &
      1.86870e-08 ,1.97458e-08 ,2.08553e-08 ,2.20175e-08 ,2.32342e-08 , &
      2.45077e-08 ,2.58401e-08 ,2.72334e-08 ,2.86900e-08 ,3.02122e-08 , &
      3.18021e-08 ,3.34624e-08 ,3.51954e-08 ,3.70037e-08 ,3.88899e-08 /)
      totplnk(101:150,14) = (/ &
      4.08568e-08 ,4.29068e-08 ,4.50429e-08 ,4.72678e-08 ,4.95847e-08 , &
      5.19963e-08 ,5.45058e-08 ,5.71161e-08 ,5.98309e-08 ,6.26529e-08 , &
      6.55857e-08 ,6.86327e-08 ,7.17971e-08 ,7.50829e-08 ,7.84933e-08 , &
      8.20323e-08 ,8.57035e-08 ,8.95105e-08 ,9.34579e-08 ,9.75488e-08 , &
      1.01788e-07 ,1.06179e-07 ,1.10727e-07 ,1.15434e-07 ,1.20307e-07 , &
      1.25350e-07 ,1.30566e-07 ,1.35961e-07 ,1.41539e-07 ,1.47304e-07 , &
      1.53263e-07 ,1.59419e-07 ,1.65778e-07 ,1.72345e-07 ,1.79124e-07 , &
      1.86122e-07 ,1.93343e-07 ,2.00792e-07 ,2.08476e-07 ,2.16400e-07 , &
      2.24568e-07 ,2.32988e-07 ,2.41666e-07 ,2.50605e-07 ,2.59813e-07 , &
      2.69297e-07 ,2.79060e-07 ,2.89111e-07 ,2.99455e-07 ,3.10099e-07 /)
      totplnk(151:181,14) = (/ &
      3.21049e-07 ,3.32311e-07 ,3.43893e-07 ,3.55801e-07 ,3.68041e-07 , &
      3.80621e-07 ,3.93547e-07 ,4.06826e-07 ,4.20465e-07 ,4.34473e-07 , &
      4.48856e-07 ,4.63620e-07 ,4.78774e-07 ,4.94325e-07 ,5.10280e-07 , &
      5.26648e-07 ,5.43436e-07 ,5.60652e-07 ,5.78302e-07 ,5.96397e-07 , &
      6.14943e-07 ,6.33949e-07 ,6.53421e-07 ,6.73370e-07 ,6.93803e-07 , &
      7.14731e-07 ,7.36157e-07 ,7.58095e-07 ,7.80549e-07 ,8.03533e-07 , &
      8.27050e-07 /)
      totplnk(1:50,15) = (/ &
      3.90483e-12 ,4.47999e-12 ,5.13122e-12 ,5.86739e-12 ,6.69829e-12 , &
      7.63467e-12 ,8.68833e-12 ,9.87221e-12 ,1.12005e-11 ,1.26885e-11 , &
      1.43534e-11 ,1.62134e-11 ,1.82888e-11 ,2.06012e-11 ,2.31745e-11 , &
      2.60343e-11 ,2.92087e-11 ,3.27277e-11 ,3.66242e-11 ,4.09334e-11 , &
      4.56935e-11 ,5.09455e-11 ,5.67338e-11 ,6.31057e-11 ,7.01127e-11 , &
      7.78096e-11 ,8.62554e-11 ,9.55130e-11 ,1.05651e-10 ,1.16740e-10 , &
      1.28858e-10 ,1.42089e-10 ,1.56519e-10 ,1.72243e-10 ,1.89361e-10 , &
      2.07978e-10 ,2.28209e-10 ,2.50173e-10 ,2.73999e-10 ,2.99820e-10 , &
      3.27782e-10 ,3.58034e-10 ,3.90739e-10 ,4.26067e-10 ,4.64196e-10 , &
      5.05317e-10 ,5.49631e-10 ,5.97347e-10 ,6.48689e-10 ,7.03891e-10 /)
      totplnk(51:100,15) = (/ &
      7.63201e-10 ,8.26876e-10 ,8.95192e-10 ,9.68430e-10 ,1.04690e-09 , &
      1.13091e-09 ,1.22079e-09 ,1.31689e-09 ,1.41957e-09 ,1.52922e-09 , &
      1.64623e-09 ,1.77101e-09 ,1.90401e-09 ,2.04567e-09 ,2.19647e-09 , &
      2.35690e-09 ,2.52749e-09 ,2.70875e-09 ,2.90127e-09 ,3.10560e-09 , &
      3.32238e-09 ,3.55222e-09 ,3.79578e-09 ,4.05375e-09 ,4.32682e-09 , &
      4.61574e-09 ,4.92128e-09 ,5.24420e-09 ,5.58536e-09 ,5.94558e-09 , &
      6.32575e-09 ,6.72678e-09 ,7.14964e-09 ,7.59526e-09 ,8.06470e-09 , &
      8.55897e-09 ,9.07916e-09 ,9.62638e-09 ,1.02018e-08 ,1.08066e-08 , &
      1.14420e-08 ,1.21092e-08 ,1.28097e-08 ,1.35446e-08 ,1.43155e-08 , &
      1.51237e-08 ,1.59708e-08 ,1.68581e-08 ,1.77873e-08 ,1.87599e-08 /)
      totplnk(101:150,15) = (/ &
      1.97777e-08 ,2.08423e-08 ,2.19555e-08 ,2.31190e-08 ,2.43348e-08 , &
      2.56045e-08 ,2.69302e-08 ,2.83140e-08 ,2.97578e-08 ,3.12636e-08 , &
      3.28337e-08 ,3.44702e-08 ,3.61755e-08 ,3.79516e-08 ,3.98012e-08 , &
      4.17265e-08 ,4.37300e-08 ,4.58143e-08 ,4.79819e-08 ,5.02355e-08 , &
      5.25777e-08 ,5.50114e-08 ,5.75393e-08 ,6.01644e-08 ,6.28896e-08 , &
      6.57177e-08 ,6.86521e-08 ,7.16959e-08 ,7.48520e-08 ,7.81239e-08 , &
      8.15148e-08 ,8.50282e-08 ,8.86675e-08 ,9.24362e-08 ,9.63380e-08 , &
      1.00376e-07 ,1.04555e-07 ,1.08878e-07 ,1.13349e-07 ,1.17972e-07 , &
      1.22751e-07 ,1.27690e-07 ,1.32793e-07 ,1.38064e-07 ,1.43508e-07 , &
      1.49129e-07 ,1.54931e-07 ,1.60920e-07 ,1.67099e-07 ,1.73473e-07 /)
      totplnk(151:181,15) = (/ &
      1.80046e-07 ,1.86825e-07 ,1.93812e-07 ,2.01014e-07 ,2.08436e-07 , &
      2.16082e-07 ,2.23957e-07 ,2.32067e-07 ,2.40418e-07 ,2.49013e-07 , &
      2.57860e-07 ,2.66963e-07 ,2.76328e-07 ,2.85961e-07 ,2.95868e-07 , &
      3.06053e-07 ,3.16524e-07 ,3.27286e-07 ,3.38345e-07 ,3.49707e-07 , &
      3.61379e-07 ,3.73367e-07 ,3.85676e-07 ,3.98315e-07 ,4.11287e-07 , &
      4.24602e-07 ,4.38265e-07 ,4.52283e-07 ,4.66662e-07 ,4.81410e-07 , &
      4.96535e-07 /)
      totplnk(1:50,16) = (/ &
      0.28639e-12 ,0.33349e-12 ,0.38764e-12 ,0.44977e-12 ,0.52093e-12 , &
      0.60231e-12 ,0.69522e-12 ,0.80111e-12 ,0.92163e-12 ,0.10586e-11 , &
      0.12139e-11 ,0.13899e-11 ,0.15890e-11 ,0.18138e-11 ,0.20674e-11 , &
      0.23531e-11 ,0.26744e-11 ,0.30352e-11 ,0.34401e-11 ,0.38936e-11 , &
      0.44011e-11 ,0.49681e-11 ,0.56010e-11 ,0.63065e-11 ,0.70919e-11 , &
      0.79654e-11 ,0.89357e-11 ,0.10012e-10 ,0.11205e-10 ,0.12526e-10 , &
      0.13986e-10 ,0.15600e-10 ,0.17380e-10 ,0.19342e-10 ,0.21503e-10 , &
      0.23881e-10 ,0.26494e-10 ,0.29362e-10 ,0.32509e-10 ,0.35958e-10 , &
      0.39733e-10 ,0.43863e-10 ,0.48376e-10 ,0.53303e-10 ,0.58679e-10 , &
      0.64539e-10 ,0.70920e-10 ,0.77864e-10 ,0.85413e-10 ,0.93615e-10 /)
      totplnk(51:100,16) = (/ &
      0.10252e-09 ,0.11217e-09 ,0.12264e-09 ,0.13397e-09 ,0.14624e-09 , &
      0.15950e-09 ,0.17383e-09 ,0.18930e-09 ,0.20599e-09 ,0.22399e-09 , &
      0.24339e-09 ,0.26427e-09 ,0.28674e-09 ,0.31090e-09 ,0.33686e-09 , &
      0.36474e-09 ,0.39466e-09 ,0.42676e-09 ,0.46115e-09 ,0.49800e-09 , &
      0.53744e-09 ,0.57964e-09 ,0.62476e-09 ,0.67298e-09 ,0.72448e-09 , &
      0.77945e-09 ,0.83809e-09 ,0.90062e-09 ,0.96725e-09 ,0.10382e-08 , &
      0.11138e-08 ,0.11941e-08 ,0.12796e-08 ,0.13704e-08 ,0.14669e-08 , &
      0.15694e-08 ,0.16781e-08 ,0.17934e-08 ,0.19157e-08 ,0.20453e-08 , &
      0.21825e-08 ,0.23278e-08 ,0.24815e-08 ,0.26442e-08 ,0.28161e-08 , &
      0.29978e-08 ,0.31898e-08 ,0.33925e-08 ,0.36064e-08 ,0.38321e-08 /)
      totplnk(101:150,16) = (/ &
      0.40700e-08 ,0.43209e-08 ,0.45852e-08 ,0.48636e-08 ,0.51567e-08 , &
      0.54652e-08 ,0.57897e-08 ,0.61310e-08 ,0.64897e-08 ,0.68667e-08 , &
      0.72626e-08 ,0.76784e-08 ,0.81148e-08 ,0.85727e-08 ,0.90530e-08 , &
      0.95566e-08 ,0.10084e-07 ,0.10638e-07 ,0.11217e-07 ,0.11824e-07 , &
      0.12458e-07 ,0.13123e-07 ,0.13818e-07 ,0.14545e-07 ,0.15305e-07 , &
      0.16099e-07 ,0.16928e-07 ,0.17795e-07 ,0.18699e-07 ,0.19643e-07 , &
      0.20629e-07 ,0.21656e-07 ,0.22728e-07 ,0.23845e-07 ,0.25010e-07 , &
      0.26223e-07 ,0.27487e-07 ,0.28804e-07 ,0.30174e-07 ,0.31600e-07 , &
      0.33084e-07 ,0.34628e-07 ,0.36233e-07 ,0.37902e-07 ,0.39637e-07 , &
      0.41440e-07 ,0.43313e-07 ,0.45259e-07 ,0.47279e-07 ,0.49376e-07 /)
      totplnk(151:181,16) = (/ &
      0.51552e-07 ,0.53810e-07 ,0.56153e-07 ,0.58583e-07 ,0.61102e-07 , &
      0.63713e-07 ,0.66420e-07 ,0.69224e-07 ,0.72129e-07 ,0.75138e-07 , &
      0.78254e-07 ,0.81479e-07 ,0.84818e-07 ,0.88272e-07 ,0.91846e-07 , &
      0.95543e-07 ,0.99366e-07 ,0.10332e-06 ,0.10740e-06 ,0.11163e-06 , &
      0.11599e-06 ,0.12050e-06 ,0.12515e-06 ,0.12996e-06 ,0.13493e-06 , &
      0.14005e-06 ,0.14534e-06 ,0.15080e-06 ,0.15643e-06 ,0.16224e-06 , &
      0.16823e-06 /)
      totplk16(1:50) = (/ &
      0.28481e-12 ,0.33159e-12 ,0.38535e-12 ,0.44701e-12 ,0.51763e-12 , &
      0.59836e-12 ,0.69049e-12 ,0.79549e-12 ,0.91493e-12 ,0.10506e-11 , &
      0.12045e-11 ,0.13788e-11 ,0.15758e-11 ,0.17984e-11 ,0.20493e-11 , &
      0.23317e-11 ,0.26494e-11 ,0.30060e-11 ,0.34060e-11 ,0.38539e-11 , &
      0.43548e-11 ,0.49144e-11 ,0.55387e-11 ,0.62344e-11 ,0.70086e-11 , &
      0.78692e-11 ,0.88248e-11 ,0.98846e-11 ,0.11059e-10 ,0.12358e-10 , &
      0.13794e-10 ,0.15379e-10 ,0.17128e-10 ,0.19055e-10 ,0.21176e-10 , &
      0.23508e-10 ,0.26070e-10 ,0.28881e-10 ,0.31963e-10 ,0.35339e-10 , &
      0.39034e-10 ,0.43073e-10 ,0.47484e-10 ,0.52299e-10 ,0.57548e-10 , &
      0.63267e-10 ,0.69491e-10 ,0.76261e-10 ,0.83616e-10 ,0.91603e-10 /)
      totplk16(51:100) = (/ &
      0.10027e-09 ,0.10966e-09 ,0.11983e-09 ,0.13084e-09 ,0.14275e-09 , &
      0.15562e-09 ,0.16951e-09 ,0.18451e-09 ,0.20068e-09 ,0.21810e-09 , &
      0.23686e-09 ,0.25704e-09 ,0.27875e-09 ,0.30207e-09 ,0.32712e-09 , &
      0.35400e-09 ,0.38282e-09 ,0.41372e-09 ,0.44681e-09 ,0.48223e-09 , &
      0.52013e-09 ,0.56064e-09 ,0.60392e-09 ,0.65015e-09 ,0.69948e-09 , &
      0.75209e-09 ,0.80818e-09 ,0.86794e-09 ,0.93157e-09 ,0.99929e-09 , &
      0.10713e-08 ,0.11479e-08 ,0.12293e-08 ,0.13157e-08 ,0.14074e-08 , &
      0.15047e-08 ,0.16079e-08 ,0.17172e-08 ,0.18330e-08 ,0.19557e-08 , &
      0.20855e-08 ,0.22228e-08 ,0.23680e-08 ,0.25214e-08 ,0.26835e-08 , &
      0.28546e-08 ,0.30352e-08 ,0.32257e-08 ,0.34266e-08 ,0.36384e-08 /)
      totplk16(101:150) = (/ &
      0.38615e-08 ,0.40965e-08 ,0.43438e-08 ,0.46041e-08 ,0.48779e-08 , &
      0.51658e-08 ,0.54683e-08 ,0.57862e-08 ,0.61200e-08 ,0.64705e-08 , &
      0.68382e-08 ,0.72240e-08 ,0.76285e-08 ,0.80526e-08 ,0.84969e-08 , &
      0.89624e-08 ,0.94498e-08 ,0.99599e-08 ,0.10494e-07 ,0.11052e-07 , &
      0.11636e-07 ,0.12246e-07 ,0.12884e-07 ,0.13551e-07 ,0.14246e-07 , &
      0.14973e-07 ,0.15731e-07 ,0.16522e-07 ,0.17347e-07 ,0.18207e-07 , &
      0.19103e-07 ,0.20037e-07 ,0.21011e-07 ,0.22024e-07 ,0.23079e-07 , &
      0.24177e-07 ,0.25320e-07 ,0.26508e-07 ,0.27744e-07 ,0.29029e-07 , &
      0.30365e-07 ,0.31753e-07 ,0.33194e-07 ,0.34691e-07 ,0.36246e-07 , &
      0.37859e-07 ,0.39533e-07 ,0.41270e-07 ,0.43071e-07 ,0.44939e-07 /)
      totplk16(151:181) = (/ &
      0.46875e-07 ,0.48882e-07 ,0.50961e-07 ,0.53115e-07 ,0.55345e-07 , &
      0.57655e-07 ,0.60046e-07 ,0.62520e-07 ,0.65080e-07 ,0.67728e-07 , &
      0.70466e-07 ,0.73298e-07 ,0.76225e-07 ,0.79251e-07 ,0.82377e-07 , &
      0.85606e-07 ,0.88942e-07 ,0.92386e-07 ,0.95942e-07 ,0.99612e-07 , &
      0.10340e-06 ,0.10731e-06 ,0.11134e-06 ,0.11550e-06 ,0.11979e-06 , &
      0.12421e-06 ,0.12876e-06 ,0.13346e-06 ,0.13830e-06 ,0.14328e-06 , &
      0.14841e-06 /)

   end subroutine lwavplank

!***************************************************************************
   subroutine lwavplankderiv
!***************************************************************************

      save
 
      totplnkderiv(1:50,  1) = (/ &
      2.22125e-08 ,2.23245e-08 ,2.24355e-08 ,2.25435e-08 ,2.26560e-08 , &
      2.27620e-08 ,2.28690e-08 ,2.29760e-08 ,2.30775e-08 ,2.31800e-08 , &
      2.32825e-08 ,2.33825e-08 ,2.34820e-08 ,2.35795e-08 ,2.36760e-08 , &
      2.37710e-08 ,2.38655e-08 ,2.39595e-08 ,2.40530e-08 ,2.41485e-08 , &
      2.42395e-08 ,2.43300e-08 ,2.44155e-08 ,2.45085e-08 ,2.45905e-08 , &
      2.46735e-08 ,2.47565e-08 ,2.48465e-08 ,2.49315e-08 ,2.50100e-08 , &
      2.50905e-08 ,2.51705e-08 ,2.52490e-08 ,2.53260e-08 ,2.54075e-08 , &
      2.54785e-08 ,2.55555e-08 ,2.56340e-08 ,2.57050e-08 ,2.57820e-08 , &
      2.58525e-08 ,2.59205e-08 ,2.59945e-08 ,2.60680e-08 ,2.61375e-08 , &
      2.61980e-08 ,2.62745e-08 ,2.63335e-08 ,2.63995e-08 ,2.64710e-08 /)
      totplnkderiv(51:100,  1) = (/ &
      2.65300e-08 ,2.66005e-08 ,2.66685e-08 ,2.67310e-08 ,2.67915e-08 , &
      2.68540e-08 ,2.69065e-08 ,2.69730e-08 ,2.70270e-08 ,2.70690e-08 , &
      2.71420e-08 ,2.71985e-08 ,2.72560e-08 ,2.73180e-08 ,2.73760e-08 , &
      2.74285e-08 ,2.74840e-08 ,2.75290e-08 ,2.75950e-08 ,2.76360e-08 , &
      2.76975e-08 ,2.77475e-08 ,2.78080e-08 ,2.78375e-08 ,2.79120e-08 , &
      2.79510e-08 ,2.79955e-08 ,2.80625e-08 ,2.80920e-08 ,2.81570e-08 , &
      2.81990e-08 ,2.82330e-08 ,2.82830e-08 ,2.83365e-08 ,2.83740e-08 , &
      2.84295e-08 ,2.84910e-08 ,2.85275e-08 ,2.85525e-08 ,2.86085e-08 , &
      2.86535e-08 ,2.86945e-08 ,2.87355e-08 ,2.87695e-08 ,2.88105e-08 , &
      2.88585e-08 ,2.88945e-08 ,2.89425e-08 ,2.89580e-08 ,2.90265e-08 /)
      totplnkderiv(101:150,  1) = (/ &
      2.90445e-08 ,2.90905e-08 ,2.91425e-08 ,2.91560e-08 ,2.91970e-08 , &
      2.91905e-08 ,2.92880e-08 ,2.92950e-08 ,2.93630e-08 ,2.93995e-08 , &
      2.94425e-08 ,2.94635e-08 ,2.94770e-08 ,2.95290e-08 ,2.95585e-08 , &
      2.95815e-08 ,2.95995e-08 ,2.96745e-08 ,2.96725e-08 ,2.97040e-08 , &
      2.97750e-08 ,2.97905e-08 ,2.98175e-08 ,2.98355e-08 ,2.98705e-08 , &
      2.99040e-08 ,2.99680e-08 ,2.99860e-08 ,3.00270e-08 ,3.00200e-08 , &
      3.00770e-08 ,3.00795e-08 ,3.01065e-08 ,3.01795e-08 ,3.01815e-08 , &
      3.02025e-08 ,3.02360e-08 ,3.02360e-08 ,3.03090e-08 ,3.03155e-08 , &
      3.03725e-08 ,3.03635e-08 ,3.04270e-08 ,3.04610e-08 ,3.04635e-08 , &
      3.04610e-08 ,3.05180e-08 ,3.05430e-08 ,3.05290e-08 ,3.05885e-08 /)
      totplnkderiv(151:181,  1) = (/ &
      3.05750e-08 ,3.05775e-08 ,3.06795e-08 ,3.07025e-08 ,3.07365e-08 , &
      3.07435e-08 ,3.07525e-08 ,3.07680e-08 ,3.08115e-08 ,3.07930e-08 , &
      3.08155e-08 ,3.08660e-08 ,3.08865e-08 ,3.08390e-08 ,3.09340e-08 , &
      3.09685e-08 ,3.09340e-08 ,3.09820e-08 ,3.10365e-08 ,3.10705e-08 , &
      3.10750e-08 ,3.10475e-08 ,3.11685e-08 ,3.11455e-08 ,3.11500e-08 , &
      3.11775e-08 ,3.11890e-08 ,3.12045e-08 ,3.12185e-08 ,3.12415e-08 , &
      3.12590e-08 /)
      totplnkderiv(1:50,  2) = (/ &
      4.91150e-08 ,4.97290e-08 ,5.03415e-08 ,5.09460e-08 ,5.15550e-08 , &
      5.21540e-08 ,5.27575e-08 ,5.33500e-08 ,5.39500e-08 ,5.45445e-08 , &
      5.51290e-08 ,5.57235e-08 ,5.62955e-08 ,5.68800e-08 ,5.74620e-08 , &
      5.80425e-08 ,5.86145e-08 ,5.91810e-08 ,5.97435e-08 ,6.03075e-08 , &
      6.08625e-08 ,6.14135e-08 ,6.19775e-08 ,6.25185e-08 ,6.30675e-08 , &
      6.36145e-08 ,6.41535e-08 ,6.46920e-08 ,6.52265e-08 ,6.57470e-08 , &
      6.62815e-08 ,6.68000e-08 ,6.73320e-08 ,6.78550e-08 ,6.83530e-08 , &
      6.88760e-08 ,6.93735e-08 ,6.98790e-08 ,7.03950e-08 ,7.08810e-08 , &
      7.13815e-08 ,7.18795e-08 ,7.23415e-08 ,7.28505e-08 ,7.33285e-08 , &
      7.38075e-08 ,7.42675e-08 ,7.47605e-08 ,7.52380e-08 ,7.57020e-08 /)
      totplnkderiv(51:100,  2) = (/ &
      7.61495e-08 ,7.65955e-08 ,7.70565e-08 ,7.75185e-08 ,7.79735e-08 , &
      7.83915e-08 ,7.88625e-08 ,7.93215e-08 ,7.97425e-08 ,8.02195e-08 , &
      8.05905e-08 ,8.10335e-08 ,8.14770e-08 ,8.19025e-08 ,8.22955e-08 , &
      8.27115e-08 ,8.31165e-08 ,8.35645e-08 ,8.39440e-08 ,8.43785e-08 , &
      8.47380e-08 ,8.51495e-08 ,8.55405e-08 ,8.59720e-08 ,8.63135e-08 , &
      8.67065e-08 ,8.70930e-08 ,8.74545e-08 ,8.78780e-08 ,8.82160e-08 , &
      8.85625e-08 ,8.89850e-08 ,8.93395e-08 ,8.97080e-08 ,9.00675e-08 , &
      9.04085e-08 ,9.07360e-08 ,9.11315e-08 ,9.13815e-08 ,9.18320e-08 , &
      9.21500e-08 ,9.24725e-08 ,9.28640e-08 ,9.31955e-08 ,9.35185e-08 , &
      9.38645e-08 ,9.41780e-08 ,9.45465e-08 ,9.48470e-08 ,9.51375e-08 /)
      totplnkderiv(101:150,  2) = (/ &
      9.55245e-08 ,9.57925e-08 ,9.61195e-08 ,9.64750e-08 ,9.68110e-08 , &
      9.71715e-08 ,9.74150e-08 ,9.77250e-08 ,9.79600e-08 ,9.82600e-08 , &
      9.85300e-08 ,9.88400e-08 ,9.91600e-08 ,9.95350e-08 ,9.97500e-08 , &
      1.00090e-07 ,1.00370e-07 ,1.00555e-07 ,1.00935e-07 ,1.01275e-07 , &
      1.01400e-07 ,1.01790e-07 ,1.01945e-07 ,1.02225e-07 ,1.02585e-07 , &
      1.02895e-07 ,1.03010e-07 ,1.03285e-07 ,1.03540e-07 ,1.03890e-07 , &
      1.04015e-07 ,1.04420e-07 ,1.04640e-07 ,1.04810e-07 ,1.05090e-07 , &
      1.05385e-07 ,1.05600e-07 ,1.05965e-07 ,1.06050e-07 ,1.06385e-07 , &
      1.06390e-07 ,1.06795e-07 ,1.06975e-07 ,1.07240e-07 ,1.07435e-07 , &
      1.07815e-07 ,1.07960e-07 ,1.08010e-07 ,1.08535e-07 ,1.08670e-07 /)
      totplnkderiv(151:181,  2) = (/ &
      1.08855e-07 ,1.09210e-07 ,1.09195e-07 ,1.09510e-07 ,1.09665e-07 , &
      1.09885e-07 ,1.10130e-07 ,1.10440e-07 ,1.10640e-07 ,1.10760e-07 , &
      1.11125e-07 ,1.11195e-07 ,1.11345e-07 ,1.11710e-07 ,1.11765e-07 , &
      1.11960e-07 ,1.12225e-07 ,1.12460e-07 ,1.12595e-07 ,1.12730e-07 , &
      1.12880e-07 ,1.13295e-07 ,1.13215e-07 ,1.13505e-07 ,1.13665e-07 , &
      1.13870e-07 ,1.14025e-07 ,1.14325e-07 ,1.14495e-07 ,1.14605e-07 , &
      1.14905e-07 /)
      totplnkderiv(1:50, 3) = (/ &
      4.27040e-08 ,4.35430e-08 ,4.43810e-08 ,4.52210e-08 ,4.60630e-08 , &
      4.69135e-08 ,4.77585e-08 ,4.86135e-08 ,4.94585e-08 ,5.03230e-08 , &
      5.11740e-08 ,5.20250e-08 ,5.28940e-08 ,5.37465e-08 ,5.46175e-08 , &
      5.54700e-08 ,5.63430e-08 ,5.72085e-08 ,5.80735e-08 ,5.89430e-08 , &
      5.98015e-08 ,6.06680e-08 ,6.15380e-08 ,6.24130e-08 ,6.32755e-08 , &
      6.41340e-08 ,6.50060e-08 ,6.58690e-08 ,6.67315e-08 ,6.76025e-08 , &
      6.84585e-08 ,6.93205e-08 ,7.01845e-08 ,7.10485e-08 ,7.19160e-08 , &
      7.27695e-08 ,7.36145e-08 ,7.44840e-08 ,7.53405e-08 ,7.61770e-08 , &
      7.70295e-08 ,7.78745e-08 ,7.87350e-08 ,7.95740e-08 ,8.04150e-08 , &
      8.12565e-08 ,8.20885e-08 ,8.29455e-08 ,8.37830e-08 ,8.46035e-08 /)
      totplnkderiv(51:100, 3) = (/ &
      8.54315e-08 ,8.62770e-08 ,8.70975e-08 ,8.79140e-08 ,8.87190e-08 , &
      8.95625e-08 ,9.03625e-08 ,9.11795e-08 ,9.19930e-08 ,9.27685e-08 , &
      9.36095e-08 ,9.43785e-08 ,9.52375e-08 ,9.59905e-08 ,9.67680e-08 , &
      9.75840e-08 ,9.83755e-08 ,9.91710e-08 ,9.99445e-08 ,1.00706e-07 , &
      1.01477e-07 ,1.02255e-07 ,1.03021e-07 ,1.03776e-07 ,1.04544e-07 , &
      1.05338e-07 ,1.06082e-07 ,1.06843e-07 ,1.07543e-07 ,1.08298e-07 , &
      1.09103e-07 ,1.09812e-07 ,1.10536e-07 ,1.11268e-07 ,1.12027e-07 , &
      1.12727e-07 ,1.13464e-07 ,1.14183e-07 ,1.15037e-07 ,1.15615e-07 , &
      1.16329e-07 ,1.17057e-07 ,1.17734e-07 ,1.18448e-07 ,1.19149e-07 , &
      1.19835e-07 ,1.20512e-07 ,1.21127e-07 ,1.21895e-07 ,1.22581e-07 /)
      totplnkderiv(101:150, 3) = (/ &
      1.23227e-07 ,1.23928e-07 ,1.24560e-07 ,1.25220e-07 ,1.25895e-07 , &
      1.26565e-07 ,1.27125e-07 ,1.27855e-07 ,1.28490e-07 ,1.29195e-07 , &
      1.29790e-07 ,1.30470e-07 ,1.31070e-07 ,1.31690e-07 ,1.32375e-07 , &
      1.32960e-07 ,1.33570e-07 ,1.34230e-07 ,1.34840e-07 ,1.35315e-07 , &
      1.35990e-07 ,1.36555e-07 ,1.37265e-07 ,1.37945e-07 ,1.38425e-07 , &
      1.38950e-07 ,1.39640e-07 ,1.40220e-07 ,1.40775e-07 ,1.41400e-07 , &
      1.42020e-07 ,1.42500e-07 ,1.43085e-07 ,1.43680e-07 ,1.44255e-07 , &
      1.44855e-07 ,1.45385e-07 ,1.45890e-07 ,1.46430e-07 ,1.46920e-07 , &
      1.47715e-07 ,1.48090e-07 ,1.48695e-07 ,1.49165e-07 ,1.49715e-07 , &
      1.50130e-07 ,1.50720e-07 ,1.51330e-07 ,1.51725e-07 ,1.52350e-07 /)
      totplnkderiv(151:181, 3) = (/ &
      1.52965e-07 ,1.53305e-07 ,1.53915e-07 ,1.54280e-07 ,1.54950e-07 , &
      1.55370e-07 ,1.55850e-07 ,1.56260e-07 ,1.56825e-07 ,1.57470e-07 , &
      1.57760e-07 ,1.58295e-07 ,1.58780e-07 ,1.59470e-07 ,1.59940e-07 , &
      1.60325e-07 ,1.60825e-07 ,1.61100e-07 ,1.61605e-07 ,1.62045e-07 , &
      1.62670e-07 ,1.63020e-07 ,1.63625e-07 ,1.63900e-07 ,1.64420e-07 , &
      1.64705e-07 ,1.65430e-07 ,1.65610e-07 ,1.66220e-07 ,1.66585e-07 , &
      1.66965e-07 /)
      totplnkderiv(1:50, 4) = (/ &
      3.32829e-08 ,3.41160e-08 ,3.49626e-08 ,3.58068e-08 ,3.66765e-08 , &
      3.75320e-08 ,3.84095e-08 ,3.92920e-08 ,4.01830e-08 ,4.10715e-08 , &
      4.19735e-08 ,4.28835e-08 ,4.37915e-08 ,4.47205e-08 ,4.56410e-08 , &
      4.65770e-08 ,4.75090e-08 ,4.84530e-08 ,4.93975e-08 ,5.03470e-08 , &
      5.13000e-08 ,5.22560e-08 ,5.32310e-08 ,5.41865e-08 ,5.51655e-08 , &
      5.61590e-08 ,5.71120e-08 ,5.81075e-08 ,5.91060e-08 ,6.00895e-08 , &
      6.10750e-08 ,6.20740e-08 ,6.30790e-08 ,6.40765e-08 ,6.50940e-08 , &
      6.60895e-08 ,6.71230e-08 ,6.81200e-08 ,6.91260e-08 ,7.01485e-08 , &
      7.11625e-08 ,7.21870e-08 ,7.32010e-08 ,7.42080e-08 ,7.52285e-08 , &
      7.62930e-08 ,7.73040e-08 ,7.83185e-08 ,7.93410e-08 ,8.03560e-08 /)
      totplnkderiv(51:100, 4) = (/ &
      8.14115e-08 ,8.24200e-08 ,8.34555e-08 ,8.45100e-08 ,8.55265e-08 , &
      8.65205e-08 ,8.75615e-08 ,8.85870e-08 ,8.96175e-08 ,9.07015e-08 , &
      9.16475e-08 ,9.27525e-08 ,9.37055e-08 ,9.47375e-08 ,9.57995e-08 , &
      9.67635e-08 ,9.77980e-08 ,9.87735e-08 ,9.98485e-08 ,1.00904e-07 , &
      1.01900e-07 ,1.02876e-07 ,1.03905e-07 ,1.04964e-07 ,1.05956e-07 , &
      1.06870e-07 ,1.07952e-07 ,1.08944e-07 ,1.10003e-07 ,1.10965e-07 , &
      1.11952e-07 ,1.12927e-07 ,1.13951e-07 ,1.14942e-07 ,1.15920e-07 , &
      1.16968e-07 ,1.17877e-07 ,1.18930e-07 ,1.19862e-07 ,1.20817e-07 , &
      1.21817e-07 ,1.22791e-07 ,1.23727e-07 ,1.24751e-07 ,1.25697e-07 , &
      1.26634e-07 ,1.27593e-07 ,1.28585e-07 ,1.29484e-07 ,1.30485e-07 /)
      totplnkderiv(101:150, 4) = (/ &
      1.31363e-07 ,1.32391e-07 ,1.33228e-07 ,1.34155e-07 ,1.35160e-07 , &
      1.36092e-07 ,1.37070e-07 ,1.37966e-07 ,1.38865e-07 ,1.39740e-07 , &
      1.40770e-07 ,1.41620e-07 ,1.42605e-07 ,1.43465e-07 ,1.44240e-07 , &
      1.45305e-07 ,1.46220e-07 ,1.47070e-07 ,1.47935e-07 ,1.48890e-07 , &
      1.49905e-07 ,1.50640e-07 ,1.51435e-07 ,1.52335e-07 ,1.53235e-07 , &
      1.54045e-07 ,1.54895e-07 ,1.55785e-07 ,1.56870e-07 ,1.57360e-07 , &
      1.58395e-07 ,1.59185e-07 ,1.60060e-07 ,1.60955e-07 ,1.61770e-07 , &
      1.62445e-07 ,1.63415e-07 ,1.64170e-07 ,1.65125e-07 ,1.65995e-07 , &
      1.66545e-07 ,1.67580e-07 ,1.68295e-07 ,1.69130e-07 ,1.69935e-07 , &
      1.70800e-07 ,1.71610e-07 ,1.72365e-07 ,1.73215e-07 ,1.73770e-07 /)
      totplnkderiv(151:181, 4) = (/ &
      1.74590e-07 ,1.75525e-07 ,1.76095e-07 ,1.77125e-07 ,1.77745e-07 , &
      1.78580e-07 ,1.79315e-07 ,1.80045e-07 ,1.80695e-07 ,1.81580e-07 , &
      1.82360e-07 ,1.83205e-07 ,1.84055e-07 ,1.84315e-07 ,1.85225e-07 , &
      1.85865e-07 ,1.86660e-07 ,1.87445e-07 ,1.88350e-07 ,1.88930e-07 , &
      1.89420e-07 ,1.90275e-07 ,1.90630e-07 ,1.91650e-07 ,1.92485e-07 , &
      1.93285e-07 ,1.93695e-07 ,1.94595e-07 ,1.94895e-07 ,1.95960e-07 , &
      1.96525e-07 /)
      totplnkderiv(1:50, 5) = (/ &
      2.41948e-08 ,2.49273e-08 ,2.56705e-08 ,2.64263e-08 ,2.71899e-08 , &
      2.79687e-08 ,2.87531e-08 ,2.95520e-08 ,3.03567e-08 ,3.11763e-08 , &
      3.20014e-08 ,3.28390e-08 ,3.36865e-08 ,3.45395e-08 ,3.54083e-08 , &
      3.62810e-08 ,3.71705e-08 ,3.80585e-08 ,3.89650e-08 ,3.98750e-08 , &
      4.07955e-08 ,4.17255e-08 ,4.26635e-08 ,4.36095e-08 ,4.45605e-08 , &
      4.55190e-08 ,4.64910e-08 ,4.74670e-08 ,4.84480e-08 ,4.94430e-08 , &
      5.04460e-08 ,5.14440e-08 ,5.24500e-08 ,5.34835e-08 ,5.44965e-08 , &
      5.55325e-08 ,5.65650e-08 ,5.76050e-08 ,5.86615e-08 ,5.97175e-08 , &
      6.07750e-08 ,6.18400e-08 ,6.29095e-08 ,6.39950e-08 ,6.50665e-08 , &
      6.61405e-08 ,6.72290e-08 ,6.82800e-08 ,6.94445e-08 ,7.05460e-08 /)
      totplnkderiv(51:100, 5) = (/ &
      7.16400e-08 ,7.27475e-08 ,7.38790e-08 ,7.49845e-08 ,7.61270e-08 , &
      7.72375e-08 ,7.83770e-08 ,7.95045e-08 ,8.06315e-08 ,8.17715e-08 , &
      8.29275e-08 ,8.40555e-08 ,8.52110e-08 ,8.63565e-08 ,8.75045e-08 , &
      8.86735e-08 ,8.98150e-08 ,9.09970e-08 ,9.21295e-08 ,9.32730e-08 , &
      9.44605e-08 ,9.56170e-08 ,9.67885e-08 ,9.79275e-08 ,9.91190e-08 , &
      1.00278e-07 ,1.01436e-07 ,1.02625e-07 ,1.03792e-07 ,1.04989e-07 , &
      1.06111e-07 ,1.07320e-07 ,1.08505e-07 ,1.09626e-07 ,1.10812e-07 , &
      1.11948e-07 ,1.13162e-07 ,1.14289e-07 ,1.15474e-07 ,1.16661e-07 , &
      1.17827e-07 ,1.19023e-07 ,1.20167e-07 ,1.21356e-07 ,1.22499e-07 , &
      1.23653e-07 ,1.24876e-07 ,1.25983e-07 ,1.27175e-07 ,1.28325e-07 /)
      totplnkderiv(101:150, 5) = (/ &
      1.29517e-07 ,1.30685e-07 ,1.31840e-07 ,1.33013e-07 ,1.34160e-07 , &
      1.35297e-07 ,1.36461e-07 ,1.37630e-07 ,1.38771e-07 ,1.39913e-07 , &
      1.41053e-07 ,1.42218e-07 ,1.43345e-07 ,1.44460e-07 ,1.45692e-07 , &
      1.46697e-07 ,1.47905e-07 ,1.49010e-07 ,1.50210e-07 ,1.51285e-07 , &
      1.52380e-07 ,1.53555e-07 ,1.54655e-07 ,1.55805e-07 ,1.56850e-07 , &
      1.58055e-07 ,1.59115e-07 ,1.60185e-07 ,1.61255e-07 ,1.62465e-07 , &
      1.63575e-07 ,1.64675e-07 ,1.65760e-07 ,1.66765e-07 ,1.67945e-07 , &
      1.69070e-07 ,1.70045e-07 ,1.71145e-07 ,1.72260e-07 ,1.73290e-07 , &
      1.74470e-07 ,1.75490e-07 ,1.76515e-07 ,1.77555e-07 ,1.78660e-07 , &
      1.79670e-07 ,1.80705e-07 ,1.81895e-07 ,1.82745e-07 ,1.83950e-07 /)
      totplnkderiv(151:181, 5) = (/ &
      1.84955e-07 ,1.85940e-07 ,1.87080e-07 ,1.88010e-07 ,1.89145e-07 , &
      1.90130e-07 ,1.91110e-07 ,1.92130e-07 ,1.93205e-07 ,1.94230e-07 , &
      1.95045e-07 ,1.96070e-07 ,1.97155e-07 ,1.98210e-07 ,1.99080e-07 , &
      2.00280e-07 ,2.01135e-07 ,2.02150e-07 ,2.03110e-07 ,2.04135e-07 , &
      2.05110e-07 ,2.06055e-07 ,2.07120e-07 ,2.08075e-07 ,2.08975e-07 , &
      2.09950e-07 ,2.10870e-07 ,2.11830e-07 ,2.12960e-07 ,2.13725e-07 , &
      2.14765e-07 /)
      totplnkderiv(1:50, 6) = (/ &
      1.36567e-08 ,1.41766e-08 ,1.47079e-08 ,1.52499e-08 ,1.58075e-08 , &
      1.63727e-08 ,1.69528e-08 ,1.75429e-08 ,1.81477e-08 ,1.87631e-08 , &
      1.93907e-08 ,2.00297e-08 ,2.06808e-08 ,2.13432e-08 ,2.20183e-08 , &
      2.27076e-08 ,2.34064e-08 ,2.41181e-08 ,2.48400e-08 ,2.55750e-08 , &
      2.63231e-08 ,2.70790e-08 ,2.78502e-08 ,2.86326e-08 ,2.94259e-08 , &
      3.02287e-08 ,3.10451e-08 ,3.18752e-08 ,3.27108e-08 ,3.35612e-08 , &
      3.44198e-08 ,3.52930e-08 ,3.61785e-08 ,3.70690e-08 ,3.79725e-08 , &
      3.88845e-08 ,3.98120e-08 ,4.07505e-08 ,4.16965e-08 ,4.26515e-08 , &
      4.36190e-08 ,4.45925e-08 ,4.55760e-08 ,4.65735e-08 ,4.75835e-08 , &
      4.85970e-08 ,4.96255e-08 ,5.06975e-08 ,5.16950e-08 ,5.27530e-08 /)
      totplnkderiv(51:100, 6) = (/ &
      5.38130e-08 ,5.48860e-08 ,5.59715e-08 ,5.70465e-08 ,5.81385e-08 , &
      5.92525e-08 ,6.03565e-08 ,6.14815e-08 ,6.26175e-08 ,6.37475e-08 , &
      6.48855e-08 ,6.60340e-08 ,6.71980e-08 ,6.83645e-08 ,6.95430e-08 , &
      7.07145e-08 ,7.19015e-08 ,7.30995e-08 ,7.43140e-08 ,7.55095e-08 , &
      7.67115e-08 ,7.79485e-08 ,7.91735e-08 ,8.03925e-08 ,8.16385e-08 , &
      8.28775e-08 ,8.41235e-08 ,8.53775e-08 ,8.66405e-08 ,8.78940e-08 , &
      8.91805e-08 ,9.04515e-08 ,9.17290e-08 ,9.30230e-08 ,9.43145e-08 , &
      9.56200e-08 ,9.69160e-08 ,9.82140e-08 ,9.95285e-08 ,1.00829e-07 , &
      1.02145e-07 ,1.03478e-07 ,1.04787e-07 ,1.06095e-07 ,1.07439e-07 , &
      1.08785e-07 ,1.10078e-07 ,1.11466e-07 ,1.12795e-07 ,1.14133e-07 /)
      totplnkderiv(101:150, 6) = (/ &
      1.15479e-07 ,1.16825e-07 ,1.18191e-07 ,1.19540e-07 ,1.20908e-07 , &
      1.22257e-07 ,1.23634e-07 ,1.24992e-07 ,1.26345e-07 ,1.27740e-07 , &
      1.29098e-07 ,1.30447e-07 ,1.31831e-07 ,1.33250e-07 ,1.34591e-07 , &
      1.36011e-07 ,1.37315e-07 ,1.38721e-07 ,1.40103e-07 ,1.41504e-07 , &
      1.42882e-07 ,1.44259e-07 ,1.45674e-07 ,1.46997e-07 ,1.48412e-07 , &
      1.49794e-07 ,1.51167e-07 ,1.52577e-07 ,1.53941e-07 ,1.55369e-07 , &
      1.56725e-07 ,1.58125e-07 ,1.59460e-07 ,1.60895e-07 ,1.62260e-07 , &
      1.63610e-07 ,1.65085e-07 ,1.66410e-07 ,1.67805e-07 ,1.69185e-07 , &
      1.70570e-07 ,1.71915e-07 ,1.73375e-07 ,1.74775e-07 ,1.76090e-07 , &
      1.77485e-07 ,1.78905e-07 ,1.80190e-07 ,1.81610e-07 ,1.82960e-07 /)
      totplnkderiv(151:181, 6) = (/ &
      1.84330e-07 ,1.85750e-07 ,1.87060e-07 ,1.88470e-07 ,1.89835e-07 , &
      1.91250e-07 ,1.92565e-07 ,1.93925e-07 ,1.95220e-07 ,1.96620e-07 , &
      1.98095e-07 ,1.99330e-07 ,2.00680e-07 ,2.02090e-07 ,2.03360e-07 , &
      2.04775e-07 ,2.06080e-07 ,2.07440e-07 ,2.08820e-07 ,2.10095e-07 , &
      2.11445e-07 ,2.12785e-07 ,2.14050e-07 ,2.15375e-07 ,2.16825e-07 , &
      2.18080e-07 ,2.19345e-07 ,2.20710e-07 ,2.21980e-07 ,2.23425e-07 , &
      2.24645e-07 /)
      totplnkderiv(1:50, 7) = (/ &
      7.22270e-09 ,7.55350e-09 ,7.89480e-09 ,8.24725e-09 ,8.60780e-09 , &
      8.98215e-09 ,9.36430e-09 ,9.76035e-09 ,1.01652e-08 ,1.05816e-08 , &
      1.10081e-08 ,1.14480e-08 ,1.18981e-08 ,1.23600e-08 ,1.28337e-08 , &
      1.33172e-08 ,1.38139e-08 ,1.43208e-08 ,1.48413e-08 ,1.53702e-08 , &
      1.59142e-08 ,1.64704e-08 ,1.70354e-08 ,1.76178e-08 ,1.82065e-08 , &
      1.88083e-08 ,1.94237e-08 ,2.00528e-08 ,2.06913e-08 ,2.13413e-08 , &
      2.20058e-08 ,2.26814e-08 ,2.33686e-08 ,2.40729e-08 ,2.47812e-08 , &
      2.55099e-08 ,2.62449e-08 ,2.69966e-08 ,2.77569e-08 ,2.85269e-08 , &
      2.93144e-08 ,3.01108e-08 ,3.09243e-08 ,3.17433e-08 ,3.25756e-08 , &
      3.34262e-08 ,3.42738e-08 ,3.51480e-08 ,3.60285e-08 ,3.69160e-08 /)
      totplnkderiv(51:100, 7) = (/ &
      3.78235e-08 ,3.87390e-08 ,3.96635e-08 ,4.06095e-08 ,4.15600e-08 , &
      4.25180e-08 ,4.34895e-08 ,4.44800e-08 ,4.54715e-08 ,4.64750e-08 , &
      4.74905e-08 ,4.85210e-08 ,4.95685e-08 ,5.06135e-08 ,5.16725e-08 , &
      5.27480e-08 ,5.38265e-08 ,5.49170e-08 ,5.60120e-08 ,5.71275e-08 , &
      5.82610e-08 ,5.93775e-08 ,6.05245e-08 ,6.17025e-08 ,6.28355e-08 , &
      6.40135e-08 ,6.52015e-08 ,6.63865e-08 ,6.75790e-08 ,6.88120e-08 , &
      7.00070e-08 ,7.12335e-08 ,7.24720e-08 ,7.37340e-08 ,7.49775e-08 , &
      7.62415e-08 ,7.75185e-08 ,7.87915e-08 ,8.00875e-08 ,8.13630e-08 , &
      8.26710e-08 ,8.39645e-08 ,8.53060e-08 ,8.66305e-08 ,8.79915e-08 , &
      8.93080e-08 ,9.06560e-08 ,9.19860e-08 ,9.33550e-08 ,9.47305e-08 /)
      totplnkderiv(101:150, 7) = (/ &
      9.61180e-08 ,9.74500e-08 ,9.88850e-08 ,1.00263e-07 ,1.01688e-07 , &
      1.03105e-07 ,1.04489e-07 ,1.05906e-07 ,1.07345e-07 ,1.08771e-07 , &
      1.10220e-07 ,1.11713e-07 ,1.13098e-07 ,1.14515e-07 ,1.16019e-07 , &
      1.17479e-07 ,1.18969e-07 ,1.20412e-07 ,1.21852e-07 ,1.23387e-07 , &
      1.24851e-07 ,1.26319e-07 ,1.27811e-07 ,1.29396e-07 ,1.30901e-07 , &
      1.32358e-07 ,1.33900e-07 ,1.35405e-07 ,1.36931e-07 ,1.38443e-07 , &
      1.39985e-07 ,1.41481e-07 ,1.43072e-07 ,1.44587e-07 ,1.46133e-07 , &
      1.47698e-07 ,1.49203e-07 ,1.50712e-07 ,1.52363e-07 ,1.53795e-07 , &
      1.55383e-07 ,1.56961e-07 ,1.58498e-07 ,1.60117e-07 ,1.61745e-07 , &
      1.63190e-07 ,1.64790e-07 ,1.66370e-07 ,1.67975e-07 ,1.69555e-07 /)
      totplnkderiv(151:181, 7) = (/ &
      1.71060e-07 ,1.72635e-07 ,1.74345e-07 ,1.75925e-07 ,1.77395e-07 , &
      1.78960e-07 ,1.80620e-07 ,1.82180e-07 ,1.83840e-07 ,1.85340e-07 , &
      1.86940e-07 ,1.88550e-07 ,1.90095e-07 ,1.91670e-07 ,1.93385e-07 , &
      1.94895e-07 ,1.96500e-07 ,1.98090e-07 ,1.99585e-07 ,2.01280e-07 , &
      2.02950e-07 ,2.04455e-07 ,2.06075e-07 ,2.07635e-07 ,2.09095e-07 , &
      2.10865e-07 ,2.12575e-07 ,2.14050e-07 ,2.15630e-07 ,2.17060e-07 , &
      2.18715e-07 /)
      totplnkderiv(1:50, 8) = (/ &
      4.26397e-09 ,4.48470e-09 ,4.71299e-09 ,4.94968e-09 ,5.19542e-09 , &
      5.44847e-09 ,5.71195e-09 ,5.98305e-09 ,6.26215e-09 ,6.55290e-09 , &
      6.85190e-09 ,7.15950e-09 ,7.47745e-09 ,7.80525e-09 ,8.14190e-09 , &
      8.48915e-09 ,8.84680e-09 ,9.21305e-09 ,9.59105e-09 ,9.98130e-09 , &
      1.03781e-08 ,1.07863e-08 ,1.12094e-08 ,1.16371e-08 ,1.20802e-08 , &
      1.25327e-08 ,1.29958e-08 ,1.34709e-08 ,1.39592e-08 ,1.44568e-08 , &
      1.49662e-08 ,1.54828e-08 ,1.60186e-08 ,1.65612e-08 ,1.71181e-08 , &
      1.76822e-08 ,1.82591e-08 ,1.88487e-08 ,1.94520e-08 ,2.00691e-08 , &
      2.06955e-08 ,2.13353e-08 ,2.19819e-08 ,2.26479e-08 ,2.33234e-08 , &
      2.40058e-08 ,2.47135e-08 ,2.54203e-08 ,2.61414e-08 ,2.68778e-08 /)
      totplnkderiv(51:100, 8) = (/ &
      2.76265e-08 ,2.83825e-08 ,2.91632e-08 ,2.99398e-08 ,3.07389e-08 , &
      3.15444e-08 ,3.23686e-08 ,3.31994e-08 ,3.40487e-08 ,3.49020e-08 , &
      3.57715e-08 ,3.66515e-08 ,3.75465e-08 ,3.84520e-08 ,3.93675e-08 , &
      4.02985e-08 ,4.12415e-08 ,4.21965e-08 ,4.31630e-08 ,4.41360e-08 , &
      4.51220e-08 ,4.61235e-08 ,4.71440e-08 ,4.81515e-08 ,4.91905e-08 , &
      5.02395e-08 ,5.12885e-08 ,5.23735e-08 ,5.34460e-08 ,5.45245e-08 , &
      5.56375e-08 ,5.67540e-08 ,5.78780e-08 ,5.90065e-08 ,6.01520e-08 , &
      6.13000e-08 ,6.24720e-08 ,6.36530e-08 ,6.48500e-08 ,6.60500e-08 , &
      6.72435e-08 ,6.84735e-08 ,6.97025e-08 ,7.09530e-08 ,7.21695e-08 , &
      7.34270e-08 ,7.47295e-08 ,7.59915e-08 ,7.72685e-08 ,7.85925e-08 /)
      totplnkderiv(101:150, 8) = (/ &
      7.98855e-08 ,8.12205e-08 ,8.25120e-08 ,8.38565e-08 ,8.52005e-08 , &
      8.65570e-08 ,8.79075e-08 ,8.92920e-08 ,9.06535e-08 ,9.20455e-08 , &
      9.34230e-08 ,9.48355e-08 ,9.62720e-08 ,9.76890e-08 ,9.90755e-08 , &
      1.00528e-07 ,1.01982e-07 ,1.03436e-07 ,1.04919e-07 ,1.06368e-07 , &
      1.07811e-07 ,1.09326e-07 ,1.10836e-07 ,1.12286e-07 ,1.13803e-07 , &
      1.15326e-07 ,1.16809e-07 ,1.18348e-07 ,1.19876e-07 ,1.21413e-07 , &
      1.22922e-07 ,1.24524e-07 ,1.26049e-07 ,1.27573e-07 ,1.29155e-07 , &
      1.30708e-07 ,1.32327e-07 ,1.33958e-07 ,1.35480e-07 ,1.37081e-07 , &
      1.38716e-07 ,1.40326e-07 ,1.41872e-07 ,1.43468e-07 ,1.45092e-07 , &
      1.46806e-07 ,1.48329e-07 ,1.49922e-07 ,1.51668e-07 ,1.53241e-07 /)
      totplnkderiv(151:181, 8) = (/ &
      1.54996e-07 ,1.56561e-07 ,1.58197e-07 ,1.59884e-07 ,1.61576e-07 , &
      1.63200e-07 ,1.64885e-07 ,1.66630e-07 ,1.68275e-07 ,1.69935e-07 , &
      1.71650e-07 ,1.73245e-07 ,1.75045e-07 ,1.76710e-07 ,1.78330e-07 , &
      1.79995e-07 ,1.81735e-07 ,1.83470e-07 ,1.85200e-07 ,1.86890e-07 , &
      1.88595e-07 ,1.90300e-07 ,1.91995e-07 ,1.93715e-07 ,1.95495e-07 , &
      1.97130e-07 ,1.98795e-07 ,2.00680e-07 ,2.02365e-07 ,2.04090e-07 , &
      2.05830e-07 /)
      totplnkderiv(1:50, 9) = (/ &
      1.85410e-09 ,1.96515e-09 ,2.08117e-09 ,2.20227e-09 ,2.32861e-09 , &
      2.46066e-09 ,2.59812e-09 ,2.74153e-09 ,2.89058e-09 ,3.04567e-09 , &
      3.20674e-09 ,3.37442e-09 ,3.54854e-09 ,3.72892e-09 ,3.91630e-09 , &
      4.11013e-09 ,4.31150e-09 ,4.52011e-09 ,4.73541e-09 ,4.95870e-09 , &
      5.18913e-09 ,5.42752e-09 ,5.67340e-09 ,5.92810e-09 ,6.18995e-09 , &
      6.46055e-09 ,6.73905e-09 ,7.02620e-09 ,7.32260e-09 ,7.62700e-09 , &
      7.94050e-09 ,8.26370e-09 ,8.59515e-09 ,8.93570e-09 ,9.28535e-09 , &
      9.64575e-09 ,1.00154e-08 ,1.03944e-08 ,1.07839e-08 ,1.11832e-08 , &
      1.15909e-08 ,1.20085e-08 ,1.24399e-08 ,1.28792e-08 ,1.33280e-08 , &
      1.37892e-08 ,1.42573e-08 ,1.47408e-08 ,1.52345e-08 ,1.57371e-08 /)
      totplnkderiv(51:100, 9) = (/ &
      1.62496e-08 ,1.67756e-08 ,1.73101e-08 ,1.78596e-08 ,1.84161e-08 , &
      1.89869e-08 ,1.95681e-08 ,2.01632e-08 ,2.07626e-08 ,2.13800e-08 , &
      2.20064e-08 ,2.26453e-08 ,2.32970e-08 ,2.39595e-08 ,2.46340e-08 , &
      2.53152e-08 ,2.60158e-08 ,2.67235e-08 ,2.74471e-08 ,2.81776e-08 , &
      2.89233e-08 ,2.96822e-08 ,3.04488e-08 ,3.12298e-08 ,3.20273e-08 , &
      3.28304e-08 ,3.36455e-08 ,3.44765e-08 ,3.53195e-08 ,3.61705e-08 , &
      3.70385e-08 ,3.79155e-08 ,3.88065e-08 ,3.97055e-08 ,4.06210e-08 , &
      4.15490e-08 ,4.24825e-08 ,4.34355e-08 ,4.43920e-08 ,4.53705e-08 , &
      4.63560e-08 ,4.73565e-08 ,4.83655e-08 ,4.93815e-08 ,5.04180e-08 , &
      5.14655e-08 ,5.25175e-08 ,5.35865e-08 ,5.46720e-08 ,5.57670e-08 /)
      totplnkderiv(101:150, 9) = (/ &
      5.68640e-08 ,5.79825e-08 ,5.91140e-08 ,6.02515e-08 ,6.13985e-08 , &
      6.25525e-08 ,6.37420e-08 ,6.49220e-08 ,6.61145e-08 ,6.73185e-08 , &
      6.85520e-08 ,6.97760e-08 ,7.10050e-08 ,7.22650e-08 ,7.35315e-08 , &
      7.48035e-08 ,7.60745e-08 ,7.73740e-08 ,7.86870e-08 ,7.99845e-08 , &
      8.13325e-08 ,8.26615e-08 ,8.40010e-08 ,8.53640e-08 ,8.67235e-08 , &
      8.80960e-08 ,8.95055e-08 ,9.08945e-08 ,9.23045e-08 ,9.37100e-08 , &
      9.51555e-08 ,9.65630e-08 ,9.80235e-08 ,9.94920e-08 ,1.00966e-07 , &
      1.02434e-07 ,1.03898e-07 ,1.05386e-07 ,1.06905e-07 ,1.08418e-07 , &
      1.09926e-07 ,1.11454e-07 ,1.13010e-07 ,1.14546e-07 ,1.16106e-07 , &
      1.17652e-07 ,1.19264e-07 ,1.20817e-07 ,1.22395e-07 ,1.24024e-07 /)
      totplnkderiv(151:181, 9) = (/ &
      1.25585e-07 ,1.27213e-07 ,1.28817e-07 ,1.30472e-07 ,1.32088e-07 , &
      1.33752e-07 ,1.35367e-07 ,1.37018e-07 ,1.38698e-07 ,1.40394e-07 , &
      1.42026e-07 ,1.43796e-07 ,1.45438e-07 ,1.47175e-07 ,1.48866e-07 , &
      1.50576e-07 ,1.52281e-07 ,1.54018e-07 ,1.55796e-07 ,1.57515e-07 , &
      1.59225e-07 ,1.60989e-07 ,1.62754e-07 ,1.64532e-07 ,1.66285e-07 , &
      1.68070e-07 ,1.69870e-07 ,1.71625e-07 ,1.73440e-07 ,1.75275e-07 , &
      1.77040e-07 /)
      totplnkderiv(1:50,10) = (/ &
      7.14917e-10 ,7.64833e-10 ,8.17460e-10 ,8.72980e-10 ,9.31380e-10 , &
      9.92940e-10 ,1.05746e-09 ,1.12555e-09 ,1.19684e-09 ,1.27162e-09 , &
      1.35001e-09 ,1.43229e-09 ,1.51815e-09 ,1.60831e-09 ,1.70271e-09 , &
      1.80088e-09 ,1.90365e-09 ,2.01075e-09 ,2.12261e-09 ,2.23924e-09 , &
      2.36057e-09 ,2.48681e-09 ,2.61814e-09 ,2.75506e-09 ,2.89692e-09 , &
      3.04423e-09 ,3.19758e-09 ,3.35681e-09 ,3.52113e-09 ,3.69280e-09 , &
      3.86919e-09 ,4.05205e-09 ,4.24184e-09 ,4.43877e-09 ,4.64134e-09 , &
      4.85088e-09 ,5.06670e-09 ,5.29143e-09 ,5.52205e-09 ,5.75980e-09 , &
      6.00550e-09 ,6.25840e-09 ,6.51855e-09 ,6.78800e-09 ,7.06435e-09 , &
      7.34935e-09 ,7.64220e-09 ,7.94470e-09 ,8.25340e-09 ,8.57030e-09 /)
      totplnkderiv(51:100,10) = (/ &
      8.89680e-09 ,9.23255e-09 ,9.57770e-09 ,9.93045e-09 ,1.02932e-08 , &
      1.06649e-08 ,1.10443e-08 ,1.14348e-08 ,1.18350e-08 ,1.22463e-08 , &
      1.26679e-08 ,1.30949e-08 ,1.35358e-08 ,1.39824e-08 ,1.44425e-08 , &
      1.49126e-08 ,1.53884e-08 ,1.58826e-08 ,1.63808e-08 ,1.68974e-08 , &
      1.74159e-08 ,1.79447e-08 ,1.84886e-08 ,1.90456e-08 ,1.96124e-08 , &
      2.01863e-08 ,2.07737e-08 ,2.13720e-08 ,2.19837e-08 ,2.26044e-08 , &
      2.32396e-08 ,2.38856e-08 ,2.45344e-08 ,2.52055e-08 ,2.58791e-08 , &
      2.65706e-08 ,2.72758e-08 ,2.79852e-08 ,2.87201e-08 ,2.94518e-08 , &
      3.02063e-08 ,3.09651e-08 ,3.17357e-08 ,3.25235e-08 ,3.33215e-08 , &
      3.41285e-08 ,3.49485e-08 ,3.57925e-08 ,3.66330e-08 ,3.74765e-08 /)
      totplnkderiv(101:150,10) = (/ &
      3.83675e-08 ,3.92390e-08 ,4.01330e-08 ,4.10340e-08 ,4.19585e-08 , &
      4.28815e-08 ,4.38210e-08 ,4.47770e-08 ,4.57575e-08 ,4.67325e-08 , &
      4.77170e-08 ,4.87205e-08 ,4.97410e-08 ,5.07620e-08 ,5.18180e-08 , &
      5.28540e-08 ,5.39260e-08 ,5.50035e-08 ,5.60885e-08 ,5.71900e-08 , &
      5.82940e-08 ,5.94380e-08 ,6.05690e-08 ,6.17185e-08 ,6.28860e-08 , &
      6.40670e-08 ,6.52300e-08 ,6.64225e-08 ,6.76485e-08 ,6.88715e-08 , &
      7.00750e-08 ,7.13760e-08 ,7.25910e-08 ,7.38860e-08 ,7.51290e-08 , &
      7.64420e-08 ,7.77550e-08 ,7.90725e-08 ,8.03825e-08 ,8.17330e-08 , &
      8.30810e-08 ,8.44330e-08 ,8.57720e-08 ,8.72115e-08 ,8.85800e-08 , &
      8.99945e-08 ,9.13905e-08 ,9.28345e-08 ,9.42665e-08 ,9.56765e-08 /)
      totplnkderiv(151:181,10) = (/ &
      9.72000e-08 ,9.86780e-08 ,1.00105e-07 ,1.01616e-07 ,1.03078e-07 , &
      1.04610e-07 ,1.06154e-07 ,1.07639e-07 ,1.09242e-07 ,1.10804e-07 , &
      1.12384e-07 ,1.13871e-07 ,1.15478e-07 ,1.17066e-07 ,1.18703e-07 , &
      1.20294e-07 ,1.21930e-07 ,1.23543e-07 ,1.25169e-07 ,1.26806e-07 , &
      1.28503e-07 ,1.30233e-07 ,1.31834e-07 ,1.33596e-07 ,1.35283e-07 , &
      1.36947e-07 ,1.38594e-07 ,1.40362e-07 ,1.42131e-07 ,1.43823e-07 , &
      1.45592e-07 /)
      totplnkderiv(1:50,11) = (/ &
      2.25919e-10 ,2.43810e-10 ,2.62866e-10 ,2.83125e-10 ,3.04676e-10 , &
      3.27536e-10 ,3.51796e-10 ,3.77498e-10 ,4.04714e-10 ,4.33528e-10 , &
      4.64000e-10 ,4.96185e-10 ,5.30165e-10 ,5.65999e-10 ,6.03749e-10 , &
      6.43579e-10 ,6.85479e-10 ,7.29517e-10 ,7.75810e-10 ,8.24440e-10 , &
      8.75520e-10 ,9.29065e-10 ,9.85175e-10 ,1.04405e-09 ,1.10562e-09 , &
      1.17005e-09 ,1.23742e-09 ,1.30780e-09 ,1.38141e-09 ,1.45809e-09 , &
      1.53825e-09 ,1.62177e-09 ,1.70884e-09 ,1.79942e-09 ,1.89390e-09 , &
      1.99205e-09 ,2.09429e-09 ,2.20030e-09 ,2.31077e-09 ,2.42510e-09 , &
      2.54410e-09 ,2.66754e-09 ,2.79529e-09 ,2.92777e-09 ,3.06498e-09 , &
      3.20691e-09 ,3.35450e-09 ,3.50653e-09 ,3.66427e-09 ,3.82723e-09 /)
      totplnkderiv(51:100,11) = (/ &
      3.99549e-09 ,4.16911e-09 ,4.34892e-09 ,4.53415e-09 ,4.72504e-09 , &
      4.92197e-09 ,5.12525e-09 ,5.33485e-09 ,5.55085e-09 ,5.77275e-09 , &
      6.00105e-09 ,6.23650e-09 ,6.47855e-09 ,6.72735e-09 ,6.98325e-09 , &
      7.24695e-09 ,7.51730e-09 ,7.79480e-09 ,8.07975e-09 ,8.37170e-09 , &
      8.67195e-09 ,8.98050e-09 ,9.29575e-09 ,9.61950e-09 ,9.95150e-09 , &
      1.02912e-08 ,1.06397e-08 ,1.09964e-08 ,1.13611e-08 ,1.17348e-08 , &
      1.21158e-08 ,1.25072e-08 ,1.29079e-08 ,1.33159e-08 ,1.37342e-08 , &
      1.41599e-08 ,1.45966e-08 ,1.50438e-08 ,1.54964e-08 ,1.59605e-08 , &
      1.64337e-08 ,1.69189e-08 ,1.74134e-08 ,1.79136e-08 ,1.84272e-08 , &
      1.89502e-08 ,1.94845e-08 ,2.00248e-08 ,2.05788e-08 ,2.11455e-08 /)
      totplnkderiv(101:150,11) = (/ &
      2.17159e-08 ,2.23036e-08 ,2.28983e-08 ,2.35033e-08 ,2.41204e-08 , &
      2.47485e-08 ,2.53860e-08 ,2.60331e-08 ,2.66891e-08 ,2.73644e-08 , &
      2.80440e-08 ,2.87361e-08 ,2.94412e-08 ,3.01560e-08 ,3.08805e-08 , &
      3.16195e-08 ,3.23690e-08 ,3.31285e-08 ,3.39015e-08 ,3.46820e-08 , &
      3.54770e-08 ,3.62805e-08 ,3.70960e-08 ,3.79295e-08 ,3.87715e-08 , &
      3.96185e-08 ,4.04860e-08 ,4.13600e-08 ,4.22500e-08 ,4.31490e-08 , &
      4.40610e-08 ,4.49810e-08 ,4.59205e-08 ,4.68650e-08 ,4.78260e-08 , &
      4.87970e-08 ,4.97790e-08 ,5.07645e-08 ,5.17730e-08 ,5.27960e-08 , &
      5.38285e-08 ,5.48650e-08 ,5.59205e-08 ,5.69960e-08 ,5.80690e-08 , &
      5.91570e-08 ,6.02640e-08 ,6.13750e-08 ,6.25015e-08 ,6.36475e-08 /)
      totplnkderiv(151:181,11) = (/ &
      6.47950e-08 ,6.59510e-08 ,6.71345e-08 ,6.83175e-08 ,6.95250e-08 , &
      7.07325e-08 ,7.19490e-08 ,7.31880e-08 ,7.44315e-08 ,7.56880e-08 , &
      7.69500e-08 ,7.82495e-08 ,7.95330e-08 ,8.08450e-08 ,8.21535e-08 , &
      8.34860e-08 ,8.48330e-08 ,8.61795e-08 ,8.75480e-08 ,8.89235e-08 , &
      9.03060e-08 ,9.17045e-08 ,9.31140e-08 ,9.45240e-08 ,9.59720e-08 , &
      9.74140e-08 ,9.88825e-08 ,1.00347e-07 ,1.01825e-07 ,1.03305e-07 , &
      1.04826e-07 /)
      totplnkderiv(1:50,12) = (/ &
      2.91689e-11 ,3.20300e-11 ,3.51272e-11 ,3.84803e-11 ,4.21014e-11 , &
      4.60107e-11 ,5.02265e-11 ,5.47685e-11 ,5.96564e-11 ,6.49111e-11 , &
      7.05522e-11 ,7.66060e-11 ,8.30974e-11 ,9.00441e-11 ,9.74820e-11 , &
      1.05435e-10 ,1.13925e-10 ,1.22981e-10 ,1.32640e-10 ,1.42933e-10 , &
      1.53882e-10 ,1.65527e-10 ,1.77903e-10 ,1.91054e-10 ,2.05001e-10 , &
      2.19779e-10 ,2.35448e-10 ,2.52042e-10 ,2.69565e-10 ,2.88128e-10 , &
      3.07714e-10 ,3.28370e-10 ,3.50238e-10 ,3.73235e-10 ,3.97433e-10 , &
      4.22964e-10 ,4.49822e-10 ,4.78042e-10 ,5.07721e-10 ,5.38915e-10 , &
      5.71610e-10 ,6.05916e-10 ,6.41896e-10 ,6.79600e-10 ,7.19110e-10 , &
      7.60455e-10 ,8.03625e-10 ,8.48870e-10 ,8.96080e-10 ,9.45490e-10 /)
      totplnkderiv(51:100,12) = (/ &
      9.96930e-10 ,1.05071e-09 ,1.10679e-09 ,1.16521e-09 ,1.22617e-09 , &
      1.28945e-09 ,1.35554e-09 ,1.42427e-09 ,1.49574e-09 ,1.56984e-09 , &
      1.64695e-09 ,1.72715e-09 ,1.81034e-09 ,1.89656e-09 ,1.98613e-09 , &
      2.07898e-09 ,2.17515e-09 ,2.27498e-09 ,2.37826e-09 ,2.48517e-09 , &
      2.59566e-09 ,2.71004e-09 ,2.82834e-09 ,2.95078e-09 ,3.07686e-09 , &
      3.20739e-09 ,3.34232e-09 ,3.48162e-09 ,3.62515e-09 ,3.77337e-09 , &
      3.92614e-09 ,4.08317e-09 ,4.24567e-09 ,4.41272e-09 ,4.58524e-09 , &
      4.76245e-09 ,4.94450e-09 ,5.13235e-09 ,5.32535e-09 ,5.52415e-09 , &
      5.72770e-09 ,5.93815e-09 ,6.15315e-09 ,6.37525e-09 ,6.60175e-09 , &
      6.83485e-09 ,7.07490e-09 ,7.32060e-09 ,7.57225e-09 ,7.83035e-09 /)
      totplnkderiv(101:150,12) = (/ &
      8.09580e-09 ,8.36620e-09 ,8.64410e-09 ,8.93110e-09 ,9.22170e-09 , &
      9.52055e-09 ,9.82595e-09 ,1.01399e-08 ,1.04613e-08 ,1.07878e-08 , &
      1.11223e-08 ,1.14667e-08 ,1.18152e-08 ,1.21748e-08 ,1.25410e-08 , &
      1.29147e-08 ,1.32948e-08 ,1.36858e-08 ,1.40827e-08 ,1.44908e-08 , &
      1.49040e-08 ,1.53284e-08 ,1.57610e-08 ,1.61995e-08 ,1.66483e-08 , &
      1.71068e-08 ,1.75714e-08 ,1.80464e-08 ,1.85337e-08 ,1.90249e-08 , &
      1.95309e-08 ,2.00407e-08 ,2.05333e-08 ,2.10929e-08 ,2.16346e-08 , &
      2.21829e-08 ,2.27402e-08 ,2.33112e-08 ,2.38922e-08 ,2.44802e-08 , &
      2.50762e-08 ,2.56896e-08 ,2.63057e-08 ,2.69318e-08 ,2.75705e-08 , &
      2.82216e-08 ,2.88787e-08 ,2.95505e-08 ,3.02335e-08 ,3.09215e-08 /)
      totplnkderiv(151:181,12) = (/ &
      3.16235e-08 ,3.23350e-08 ,3.30590e-08 ,3.37960e-08 ,3.45395e-08 , &
      3.52955e-08 ,3.60615e-08 ,3.68350e-08 ,3.76265e-08 ,3.84255e-08 , &
      3.92400e-08 ,4.00485e-08 ,4.08940e-08 ,4.17310e-08 ,4.25860e-08 , &
      4.34585e-08 ,4.43270e-08 ,4.52220e-08 ,4.61225e-08 ,4.70345e-08 , &
      4.79560e-08 ,4.89000e-08 ,4.98445e-08 ,5.07985e-08 ,5.17705e-08 , &
      5.27575e-08 ,5.37420e-08 ,5.47495e-08 ,5.57725e-08 ,5.68105e-08 , &
      5.78395e-08 /)
      totplnkderiv(1:50,13) = (/ &
      5.47482e-12 ,6.09637e-12 ,6.77874e-12 ,7.52703e-12 ,8.34784e-12 , &
      9.24486e-12 ,1.02246e-11 ,1.12956e-11 ,1.24615e-11 ,1.37321e-11 , &
      1.51131e-11 ,1.66129e-11 ,1.82416e-11 ,2.00072e-11 ,2.19187e-11 , &
      2.39828e-11 ,2.62171e-11 ,2.86290e-11 ,3.12283e-11 ,3.40276e-11 , &
      3.70433e-11 ,4.02847e-11 ,4.37738e-11 ,4.75070e-11 ,5.15119e-11 , &
      5.58120e-11 ,6.04059e-11 ,6.53208e-11 ,7.05774e-11 ,7.61935e-11 , &
      8.21832e-11 ,8.85570e-11 ,9.53575e-11 ,1.02592e-10 ,1.10298e-10 , &
      1.18470e-10 ,1.27161e-10 ,1.36381e-10 ,1.46161e-10 ,1.56529e-10 , &
      1.67521e-10 ,1.79142e-10 ,1.91423e-10 ,2.04405e-10 ,2.18123e-10 , &
      2.32608e-10 ,2.47889e-10 ,2.63994e-10 ,2.80978e-10 ,2.98843e-10 /)
      totplnkderiv(51:100,13) = (/ &
      3.17659e-10 ,3.37423e-10 ,3.58206e-10 ,3.80090e-10 ,4.02996e-10 , &
      4.27065e-10 ,4.52298e-10 ,4.78781e-10 ,5.06493e-10 ,5.35576e-10 , &
      5.65942e-10 ,5.97761e-10 ,6.31007e-10 ,6.65740e-10 ,7.02095e-10 , &
      7.39945e-10 ,7.79575e-10 ,8.20845e-10 ,8.63870e-10 ,9.08680e-10 , &
      9.55385e-10 ,1.00416e-09 ,1.05464e-09 ,1.10737e-09 ,1.16225e-09 , &
      1.21918e-09 ,1.27827e-09 ,1.33988e-09 ,1.40370e-09 ,1.46994e-09 , &
      1.53850e-09 ,1.60993e-09 ,1.68382e-09 ,1.76039e-09 ,1.83997e-09 , &
      1.92182e-09 ,2.00686e-09 ,2.09511e-09 ,2.18620e-09 ,2.28034e-09 , &
      2.37753e-09 ,2.47805e-09 ,2.58193e-09 ,2.68935e-09 ,2.80064e-09 , &
      2.91493e-09 ,3.03271e-09 ,3.15474e-09 ,3.27987e-09 ,3.40936e-09 /)
      totplnkderiv(101:150,13) = (/ &
      3.54277e-09 ,3.68019e-09 ,3.82173e-09 ,3.96703e-09 ,4.11746e-09 , &
      4.27104e-09 ,4.43020e-09 ,4.59395e-09 ,4.76060e-09 ,4.93430e-09 , &
      5.11085e-09 ,5.29280e-09 ,5.48055e-09 ,5.67300e-09 ,5.86950e-09 , &
      6.07160e-09 ,6.28015e-09 ,6.49295e-09 ,6.71195e-09 ,6.93455e-09 , &
      7.16470e-09 ,7.39985e-09 ,7.64120e-09 ,7.88885e-09 ,8.13910e-09 , &
      8.39930e-09 ,8.66535e-09 ,8.93600e-09 ,9.21445e-09 ,9.49865e-09 , &
      9.78845e-09 ,1.00856e-08 ,1.04361e-08 ,1.07018e-08 ,1.10164e-08 , &
      1.13438e-08 ,1.16748e-08 ,1.20133e-08 ,1.23575e-08 ,1.27117e-08 , &
      1.30708e-08 ,1.34383e-08 ,1.38138e-08 ,1.41985e-08 ,1.45859e-08 , &
      1.49846e-08 ,1.53879e-08 ,1.58042e-08 ,1.62239e-08 ,1.66529e-08 /)
      totplnkderiv(151:181,13) = (/ &
      1.70954e-08 ,1.75422e-08 ,1.79943e-08 ,1.84537e-08 ,1.89280e-08 , &
      1.94078e-08 ,1.98997e-08 ,2.03948e-08 ,2.08956e-08 ,2.14169e-08 , &
      2.19330e-08 ,2.24773e-08 ,2.30085e-08 ,2.35676e-08 ,2.41237e-08 , &
      2.46919e-08 ,2.52720e-08 ,2.58575e-08 ,2.64578e-08 ,2.70675e-08 , &
      2.76878e-08 ,2.83034e-08 ,2.89430e-08 ,2.95980e-08 ,3.02480e-08 , &
      3.09105e-08 ,3.15980e-08 ,3.22865e-08 ,3.29755e-08 ,3.36775e-08 , &
      3.43990e-08 /)
      totplnkderiv(1:50,14) = (/ &
      1.81489e-12 ,2.03846e-12 ,2.28659e-12 ,2.56071e-12 ,2.86352e-12 , &
      3.19789e-12 ,3.56668e-12 ,3.97211e-12 ,4.41711e-12 ,4.90616e-12 , &
      5.44153e-12 ,6.02790e-12 ,6.67001e-12 ,7.37018e-12 ,8.13433e-12 , &
      8.96872e-12 ,9.87526e-12 ,1.08601e-11 ,1.19328e-11 ,1.30938e-11 , &
      1.43548e-11 ,1.57182e-11 ,1.71916e-11 ,1.87875e-11 ,2.05091e-11 , &
      2.23652e-11 ,2.43627e-11 ,2.65190e-11 ,2.88354e-11 ,3.13224e-11 , &
      3.39926e-11 ,3.68664e-11 ,3.99372e-11 ,4.32309e-11 ,4.67496e-11 , &
      5.05182e-11 ,5.45350e-11 ,5.88268e-11 ,6.34126e-11 ,6.82878e-11 , &
      7.34973e-11 ,7.90201e-11 ,8.49075e-11 ,9.11725e-11 ,9.78235e-11 , &
      1.04856e-10 ,1.12342e-10 ,1.20278e-10 ,1.28680e-10 ,1.37560e-10 /)
      totplnkderiv(51:100,14) = (/ &
      1.46953e-10 ,1.56900e-10 ,1.67401e-10 ,1.78498e-10 ,1.90161e-10 , &
      2.02523e-10 ,2.15535e-10 ,2.29239e-10 ,2.43665e-10 ,2.58799e-10 , &
      2.74767e-10 ,2.91522e-10 ,3.09141e-10 ,3.27625e-10 ,3.47011e-10 , &
      3.67419e-10 ,3.88720e-10 ,4.11066e-10 ,4.34522e-10 ,4.59002e-10 , &
      4.84657e-10 ,5.11391e-10 ,5.39524e-10 ,5.68709e-10 ,5.99240e-10 , &
      6.31295e-10 ,6.64520e-10 ,6.99200e-10 ,7.35525e-10 ,7.73135e-10 , &
      8.12440e-10 ,8.53275e-10 ,8.95930e-10 ,9.40165e-10 ,9.86260e-10 , &
      1.03423e-09 ,1.08385e-09 ,1.13567e-09 ,1.18916e-09 ,1.24469e-09 , &
      1.30262e-09 ,1.36268e-09 ,1.42479e-09 ,1.48904e-09 ,1.55557e-09 , &
      1.62478e-09 ,1.69642e-09 ,1.77023e-09 ,1.84696e-09 ,1.92646e-09 /)
      totplnkderiv(101:150,14) = (/ &
      2.00831e-09 ,2.09299e-09 ,2.18007e-09 ,2.27093e-09 ,2.36398e-09 , &
      2.46020e-09 ,2.55985e-09 ,2.66230e-09 ,2.76795e-09 ,2.87667e-09 , &
      2.98971e-09 ,3.10539e-09 ,3.22462e-09 ,3.34779e-09 ,3.47403e-09 , &
      3.60419e-09 ,3.73905e-09 ,3.87658e-09 ,4.01844e-09 ,4.16535e-09 , &
      4.31470e-09 ,4.46880e-09 ,4.62765e-09 ,4.78970e-09 ,4.95735e-09 , &
      5.12890e-09 ,5.30430e-09 ,5.48595e-09 ,5.67010e-09 ,5.86145e-09 , &
      6.05740e-09 ,6.25725e-09 ,6.46205e-09 ,6.67130e-09 ,6.88885e-09 , &
      7.10845e-09 ,7.33450e-09 ,7.56700e-09 ,7.80440e-09 ,8.04465e-09 , &
      8.29340e-09 ,8.54820e-09 ,8.80790e-09 ,9.07195e-09 ,9.34605e-09 , &
      9.62005e-09 ,9.90685e-09 ,1.01939e-08 ,1.04938e-08 ,1.07957e-08 /)
      totplnkderiv(151:181,14) = (/ &
      1.11059e-08 ,1.14208e-08 ,1.17447e-08 ,1.20717e-08 ,1.24088e-08 , &
      1.27490e-08 ,1.31020e-08 ,1.34601e-08 ,1.38231e-08 ,1.41966e-08 , &
      1.45767e-08 ,1.49570e-08 ,1.53503e-08 ,1.57496e-08 ,1.61663e-08 , &
      1.65784e-08 ,1.70027e-08 ,1.74290e-08 ,1.78730e-08 ,1.83235e-08 , &
      1.87810e-08 ,1.92418e-08 ,1.97121e-08 ,2.01899e-08 ,2.05787e-08 , &
      2.11784e-08 ,2.16824e-08 ,2.21931e-08 ,2.27235e-08 ,2.32526e-08 , &
      2.37850e-08 /)
      totplnkderiv(1:50,15) = (/ &
      5.39905e-13 ,6.11835e-13 ,6.92224e-13 ,7.81886e-13 ,8.81851e-13 , &
      9.93072e-13 ,1.11659e-12 ,1.25364e-12 ,1.40562e-12 ,1.57359e-12 , &
      1.75937e-12 ,1.96449e-12 ,2.19026e-12 ,2.43892e-12 ,2.71249e-12 , &
      3.01233e-12 ,3.34163e-12 ,3.70251e-12 ,4.09728e-12 ,4.52885e-12 , &
      4.99939e-12 ,5.51242e-12 ,6.07256e-12 ,6.68167e-12 ,7.34274e-12 , &
      8.06178e-12 ,8.84185e-12 ,9.68684e-12 ,1.06020e-11 ,1.15909e-11 , &
      1.26610e-11 ,1.38158e-11 ,1.50620e-11 ,1.64047e-11 ,1.78508e-11 , &
      1.94055e-11 ,2.10805e-11 ,2.28753e-11 ,2.48000e-11 ,2.68699e-11 , &
      2.90824e-11 ,3.14526e-11 ,3.39882e-11 ,3.67020e-11 ,3.95914e-11 , &
      4.26870e-11 ,4.59824e-11 ,4.94926e-11 ,5.32302e-11 ,5.72117e-11 /)
      totplnkderiv(51:100,15) = (/ &
      6.14475e-11 ,6.59483e-11 ,7.07393e-11 ,7.57999e-11 ,8.11980e-11 , &
      8.68920e-11 ,9.29390e-11 ,9.93335e-11 ,1.06101e-10 ,1.13263e-10 , &
      1.20827e-10 ,1.28819e-10 ,1.37255e-10 ,1.46163e-10 ,1.55547e-10 , &
      1.65428e-10 ,1.75837e-10 ,1.86816e-10 ,1.98337e-10 ,2.10476e-10 , &
      2.23218e-10 ,2.36600e-10 ,2.50651e-10 ,2.65425e-10 ,2.80895e-10 , &
      2.97102e-10 ,3.14100e-10 ,3.31919e-10 ,3.50568e-10 ,3.70064e-10 , &
      3.90464e-10 ,4.11813e-10 ,4.34111e-10 ,4.57421e-10 ,4.81717e-10 , &
      5.07039e-10 ,5.33569e-10 ,5.61137e-10 ,5.89975e-10 ,6.19980e-10 , &
      6.51170e-10 ,6.83650e-10 ,7.17520e-10 ,7.52735e-10 ,7.89390e-10 , &
      8.27355e-10 ,8.66945e-10 ,9.08020e-10 ,9.50665e-10 ,9.95055e-10 /)
      totplnkderiv(101:150,15) = (/ &
      1.04101e-09 ,1.08864e-09 ,1.13823e-09 ,1.18923e-09 ,1.24257e-09 , &
      1.29741e-09 ,1.35442e-09 ,1.41347e-09 ,1.47447e-09 ,1.53767e-09 , &
      1.60322e-09 ,1.67063e-09 ,1.74033e-09 ,1.81256e-09 ,1.88704e-09 , &
      1.96404e-09 ,2.04329e-09 ,2.12531e-09 ,2.21032e-09 ,2.29757e-09 , &
      2.38739e-09 ,2.48075e-09 ,2.57628e-09 ,2.67481e-09 ,2.77627e-09 , &
      2.88100e-09 ,2.98862e-09 ,3.09946e-09 ,3.21390e-09 ,3.33105e-09 , &
      3.45185e-09 ,3.57599e-09 ,3.70370e-09 ,3.83512e-09 ,3.96909e-09 , &
      4.10872e-09 ,4.25070e-09 ,4.39605e-09 ,4.54670e-09 ,4.70015e-09 , &
      4.85850e-09 ,5.02050e-09 ,5.18655e-09 ,5.35815e-09 ,5.53180e-09 , &
      5.71225e-09 ,5.89495e-09 ,6.08260e-09 ,6.27485e-09 ,6.47345e-09 /)
      totplnkderiv(151:181,15) = (/ &
      6.67520e-09 ,6.88310e-09 ,7.09400e-09 ,7.31140e-09 ,7.53350e-09 , &
      7.76040e-09 ,7.99215e-09 ,8.22850e-09 ,8.47235e-09 ,8.71975e-09 , &
      8.97360e-09 ,9.23365e-09 ,9.49950e-09 ,9.76965e-09 ,1.00441e-08 , &
      1.03270e-08 ,1.06158e-08 ,1.09112e-08 ,1.12111e-08 ,1.15172e-08 , &
      1.18263e-08 ,1.21475e-08 ,1.24735e-08 ,1.28027e-08 ,1.32023e-08 , &
      1.34877e-08 ,1.38399e-08 ,1.42000e-08 ,1.45625e-08 ,1.49339e-08 , &
      1.53156e-08 /)
      totplnkderiv(1:50,16) = (/ &
      4.38799e-14 ,5.04835e-14 ,5.79773e-14 ,6.64627e-14 ,7.60706e-14 , &
      8.69213e-14 ,9.91554e-14 ,1.12932e-13 ,1.28419e-13 ,1.45809e-13 , &
      1.65298e-13 ,1.87109e-13 ,2.11503e-13 ,2.38724e-13 ,2.69058e-13 , &
      3.02878e-13 ,3.40423e-13 ,3.82128e-13 ,4.28390e-13 ,4.79625e-13 , &
      5.36292e-13 ,5.98933e-13 ,6.68066e-13 ,7.44216e-13 ,8.28159e-13 , &
      9.20431e-13 ,1.02180e-12 ,1.13307e-12 ,1.25504e-12 ,1.38863e-12 , &
      1.53481e-12 ,1.69447e-12 ,1.86896e-12 ,2.05903e-12 ,2.26637e-12 , &
      2.49193e-12 ,2.73736e-12 ,3.00416e-12 ,3.29393e-12 ,3.60781e-12 , &
      3.94805e-12 ,4.31675e-12 ,4.71543e-12 ,5.14627e-12 ,5.61226e-12 , &
      6.11456e-12 ,6.65585e-12 ,7.23969e-12 ,7.86811e-12 ,8.54456e-12 /)
      totplnkderiv(51:100,16) = (/ &
      9.27075e-12 ,1.00516e-11 ,1.08898e-11 ,1.17884e-11 ,1.27514e-11 , &
      1.37839e-11 ,1.48893e-11 ,1.60716e-11 ,1.73333e-11 ,1.86849e-11 , &
      2.01237e-11 ,2.16610e-11 ,2.33001e-11 ,2.50440e-11 ,2.69035e-11 , &
      2.88827e-11 ,3.09881e-11 ,3.32234e-11 ,3.55981e-11 ,3.81193e-11 , &
      4.07946e-11 ,4.36376e-11 ,4.66485e-11 ,4.98318e-11 ,5.32080e-11 , &
      5.67754e-11 ,6.05524e-11 ,6.45450e-11 ,6.87639e-11 ,7.32160e-11 , &
      7.79170e-11 ,8.28780e-11 ,8.81045e-11 ,9.36200e-11 ,9.94280e-11 , &
      1.05545e-10 ,1.11982e-10 ,1.18752e-10 ,1.25866e-10 ,1.33350e-10 , &
      1.41210e-10 ,1.49469e-10 ,1.58143e-10 ,1.67233e-10 ,1.76760e-10 , &
      1.86758e-10 ,1.97236e-10 ,2.08227e-10 ,2.19723e-10 ,2.31737e-10 /)
      totplnkderiv(101:150,16) = (/ &
      2.44329e-10 ,2.57503e-10 ,2.71267e-10 ,2.85647e-10 ,3.00706e-10 , &
      3.16391e-10 ,3.32807e-10 ,3.49887e-10 ,3.67748e-10 ,3.86369e-10 , &
      4.05746e-10 ,4.25984e-10 ,4.47060e-10 ,4.68993e-10 ,4.91860e-10 , &
      5.15601e-10 ,5.40365e-10 ,5.66085e-10 ,5.92855e-10 ,6.20640e-10 , &
      6.49605e-10 ,6.79585e-10 ,7.10710e-10 ,7.43145e-10 ,7.76805e-10 , &
      8.11625e-10 ,8.47800e-10 ,8.85300e-10 ,9.24220e-10 ,9.64550e-10 , &
      1.00623e-09 ,1.04957e-09 ,1.09429e-09 ,1.14079e-09 ,1.18882e-09 , &
      1.23848e-09 ,1.28986e-09 ,1.34301e-09 ,1.39796e-09 ,1.45493e-09 , &
      1.51372e-09 ,1.57440e-09 ,1.63702e-09 ,1.70173e-09 ,1.76874e-09 , &
      1.83753e-09 ,1.90898e-09 ,1.98250e-09 ,2.05836e-09 ,2.13646e-09 /)
      totplnkderiv(151:181,16) = (/ &
      2.21710e-09 ,2.30027e-09 ,2.38591e-09 ,2.47432e-09 ,2.56503e-09 , &
      2.65878e-09 ,2.75516e-09 ,2.85432e-09 ,2.95688e-09 ,3.06201e-09 , &
      3.17023e-09 ,3.28153e-09 ,3.39604e-09 ,3.51391e-09 ,3.63517e-09 , &
      3.75955e-09 ,3.88756e-09 ,4.01880e-09 ,4.15405e-09 ,4.29255e-09 , &
      4.43535e-09 ,4.58145e-09 ,4.73165e-09 ,4.88560e-09 ,5.04390e-09 , &
      5.20630e-09 ,5.37255e-09 ,5.54355e-09 ,5.71915e-09 ,5.89855e-09 , &
      6.08280e-09 /)
      totplk16deriv(1:50) = (/ &
      4.35811e-14 ,5.01270e-14 ,5.75531e-14 ,6.59588e-14 ,7.54735e-14 , &
      8.62147e-14 ,9.83225e-14 ,1.11951e-13 ,1.27266e-13 ,1.44456e-13 , &
      1.63715e-13 ,1.85257e-13 ,2.09343e-13 ,2.36209e-13 ,2.66136e-13 , &
      2.99486e-13 ,3.36493e-13 ,3.77582e-13 ,4.23146e-13 ,4.73578e-13 , &
      5.29332e-13 ,5.90936e-13 ,6.58891e-13 ,7.33710e-13 ,8.16135e-13 , &
      9.06705e-13 ,1.00614e-12 ,1.11524e-12 ,1.23477e-12 ,1.36561e-12 , &
      1.50871e-12 ,1.66488e-12 ,1.83552e-12 ,2.02123e-12 ,2.22375e-12 , &
      2.44389e-12 ,2.68329e-12 ,2.94338e-12 ,3.22570e-12 ,3.53129e-12 , &
      3.86236e-12 ,4.22086e-12 ,4.60827e-12 ,5.02666e-12 ,5.47890e-12 , &
      5.96595e-12 ,6.49057e-12 ,7.05592e-12 ,7.66401e-12 ,8.31821e-12 /)
      totplk16deriv(51:100) = (/ &
      9.01998e-12 ,9.77390e-12 ,1.05826e-11 ,1.14491e-11 ,1.23769e-11 , &
      1.33709e-11 ,1.44341e-11 ,1.55706e-11 ,1.67821e-11 ,1.80793e-11 , &
      1.94586e-11 ,2.09316e-11 ,2.25007e-11 ,2.41685e-11 ,2.59454e-11 , &
      2.78356e-11 ,2.98440e-11 ,3.19744e-11 ,3.42355e-11 ,3.66340e-11 , &
      3.91772e-11 ,4.18773e-11 ,4.47339e-11 ,4.77509e-11 ,5.09490e-11 , &
      5.43240e-11 ,5.78943e-11 ,6.16648e-11 ,6.56445e-11 ,6.98412e-11 , &
      7.42680e-11 ,7.89335e-11 ,8.38450e-11 ,8.90220e-11 ,9.44695e-11 , &
      1.00197e-10 ,1.06221e-10 ,1.12550e-10 ,1.19193e-10 ,1.26175e-10 , &
      1.33498e-10 ,1.41188e-10 ,1.49251e-10 ,1.57693e-10 ,1.66530e-10 , &
      1.75798e-10 ,1.85495e-10 ,1.95661e-10 ,2.06275e-10 ,2.17357e-10 /)
      totplk16deriv(101:150) = (/ &
      2.28959e-10 ,2.41085e-10 ,2.53739e-10 ,2.66944e-10 ,2.80755e-10 , &
      2.95121e-10 ,3.10141e-10 ,3.25748e-10 ,3.42057e-10 ,3.59026e-10 , &
      3.76668e-10 ,3.95066e-10 ,4.14211e-10 ,4.34111e-10 ,4.54818e-10 , &
      4.76295e-10 ,4.98681e-10 ,5.21884e-10 ,5.46000e-10 ,5.71015e-10 , &
      5.97065e-10 ,6.23965e-10 ,6.51865e-10 ,6.80905e-10 ,7.11005e-10 , &
      7.42100e-10 ,7.74350e-10 ,8.07745e-10 ,8.42355e-10 ,8.78185e-10 , &
      9.15130e-10 ,9.53520e-10 ,9.93075e-10 ,1.03415e-09 ,1.07649e-09 , &
      1.12021e-09 ,1.16539e-09 ,1.21207e-09 ,1.26025e-09 ,1.31014e-09 , &
      1.36156e-09 ,1.41453e-09 ,1.46909e-09 ,1.52540e-09 ,1.58368e-09 , &
      1.64334e-09 ,1.70527e-09 ,1.76888e-09 ,1.83442e-09 ,1.90182e-09 /)
      totplk16deriv(151:181) = (/ &
      1.97128e-09 ,2.04281e-09 ,2.11635e-09 ,2.19219e-09 ,2.26979e-09 , &
      2.34989e-09 ,2.43219e-09 ,2.51660e-09 ,2.60396e-09 ,2.69317e-09 , &
      2.78501e-09 ,2.87927e-09 ,2.97600e-09 ,3.07548e-09 ,3.17772e-09 , &
      3.28235e-09 ,3.38982e-09 ,3.49985e-09 ,3.61307e-09 ,3.72883e-09 , &
      3.84805e-09 ,3.96975e-09 ,4.09465e-09 ,4.22240e-09 ,4.35370e-09 , &
      4.48800e-09 ,4.62535e-09 ,4.76640e-09 ,4.91110e-09 ,5.05850e-09 , &
      5.20965e-09 /)

   end subroutine lwavplankderiv

end module rrtmg_lw_setcoef
