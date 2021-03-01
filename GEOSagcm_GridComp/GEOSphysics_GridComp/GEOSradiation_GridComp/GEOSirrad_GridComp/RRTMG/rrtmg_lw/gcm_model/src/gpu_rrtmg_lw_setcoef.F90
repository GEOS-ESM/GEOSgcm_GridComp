#include "_gpudef.inc"

     module gpu_rrtmg_lw_setcoef


      use gpu_rrtmg_lw_rtrnmc
     
      use parrrtm, only : nbndlw, mg, maxxsec, mxmol
      use rrlw_wvn, only: totplnk, totplk16, totplnkderiv, totplk16deriv
      use rrlw_vsn, only: hvrset, hnamset
      use rrlw_ref, only : chi_mlsd
     
      use gpu_rrtmg_lw_taumol
   

      implicit none

      real  _gpudev, allocatable :: taveld(:,:)           ! layer temperatures (K)
                                                      !    Dimensions: (nlayers)
      real  _gpudev, allocatable :: tzd(:,:)             ! level (interface) temperatures (K)
                                                      !    Dimensions: (0:nlayers)
      real  _gpudev, allocatable :: tboundd(:)             ! surface temperature (K)
      real  _gpudev, allocatable :: wbroadd(:,:)          ! broadening gas column density (mol/cm2)
                                                      !    Dimensions: (nlayers)
     

      real  _gpudev :: totplnkd(181,nbndlw)
      real  _gpudev :: totplk16d(181)

      real  _gpudev :: totplnkderivd(181,nbndlw)
      real  _gpudev :: totplk16derivd(181)


      contains

      ! (dmb 2012) This subroutine allocates the needed GPU arrays
      subroutine allocateGPUSetCoef( ncol, nlayers )

          integer, intent(in) :: ncol
          integer, intent(in) :: nlayers
          allocate( taveld( ncol, nlayers) )
          allocate( tzd( ncol, 0:nlayers) )
          allocate( tboundd( ncol ))
          allocate( wbroadd( ncol, nlayers) )
   
      end subroutine

      ! (dmb 2012) This subroutine deallocates the GPU arrays
      subroutine deallocateGPUSetCoef( )

          deallocate( taveld )
          deallocate( tzd )
          deallocate( tboundd)
          deallocate( wbroadd)
      
      end subroutine

      ! (dmb 2012) Copy the needed reference data from the CPU to the GPU
      subroutine copyGPUSetCoef()

        totplnkd = totplnk
        totplk16d = totplk16
        totplnkderivd = totplnkderiv
        totplk16derivd = totplk16deriv

      end subroutine

!----------------------------------------------------------------------------
      _gpuker subroutine setcoefg(ncol, nlayers, istart, idrv)
!----------------------------------------------------------------------------
!
!  Purpose:  For a given atmosphere, calculate the indices and
!  fractions related to the pressure and temperature interpolations.
!  Also calculate the values of the integrated Planck functions 
!  for each band at the level and layer temperatures.

! ------- Declarations -------

! ----- Input -----
      integer , value, intent(in) :: ncol
      integer , value, intent(in) :: nlayers         ! total number of layers
      integer , value, intent(in) :: istart          ! beginning band of calculation
      integer , value, intent(in) :: idrv            ! Planck derivative option flag

     
                                                        

! ----- Local -----
      integer  :: indbound, indlev0
      integer  :: lay, indlay, indlev, iband
      integer  :: jp1
      real  :: stpfac, tbndfrac, t0frac, tlayfrac, tlevfrac
      real  :: dbdtlev, dbdtlay
      real  :: plog, fp, ft, ft1, water, scalefac, factor, compfp
      integer  :: iplon
      real  :: wv, lcoldry

#ifdef _CUDA
      iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      if (iplon <= ncol) then
#else
      do iplon = 1, ncol
#endif

      stpfac = 296. /1013. 

      indbound = tboundd(iplon) - 159. 
      if (indbound .lt. 1) then
         indbound = 1
      elseif (indbound .gt. 180) then
         indbound = 180
      endif
      tbndfrac = tboundd(iplon) - 159.  - float(indbound)
      indlev0 = tzd(iplon, 0) - 159. 
      if (indlev0 .lt. 1) then
         indlev0 = 1
      elseif (indlev0 .gt. 180) then
         indlev0 = 180
      endif
      t0frac = tzd(iplon, 0) - 159.  - float(indlev0)
      laytrop(iplon) = 0

! Begin layer loop 
!  Calculate the integrated Planck functions for each band at the
!  surface, level, and layer temperatures.
      do lay = 1, nlayers
         indlay = taveld(iplon, lay) - 159. 
         lcoldry = coldry( iplon, lay) 
         wv = colh2o(iplon, lay) * lcoldry
         if (indlay .lt. 1) then
            indlay = 1
         elseif (indlay .gt. 180) then
            indlay = 180
         endif
         tlayfrac = taveld(iplon, lay) - 159.  - float(indlay)
         indlev = tzd(iplon, lay) - 159. 
         if (indlev .lt. 1) then
            indlev = 1
         elseif (indlev .gt. 180) then
            indlev = 180
         endif
         tlevfrac = tzd(iplon, lay) - 159.  - float(indlev)

! Begin spectral band loop 
         do iband = 1, 15
            if (lay.eq.1) then
               dbdtlev = totplnkd(indbound+1,iband) - totplnkd(indbound,iband)
               plankbndd(iplon, iband) = semissd(iplon, iband) * &
                   (totplnkd(indbound,iband) + tbndfrac * dbdtlev)
               dbdtlev = totplnkd(indlev0+1,iband)-totplnkd(indlev0,iband)
               planklevd(iplon, 0,iband) = totplnkd(indlev0,iband) + t0frac * dbdtlev
               if (idrv .eq. 1) then 
                  dbdtlev = totplnkderivd(indbound+1,iband) - totplnkderivd(indbound,iband)
                  dplankbnd_dtd(iplon, iband) = semissd(iplon, iband) * &
                      (totplnkderivd(indbound,iband) + tbndfrac * dbdtlev)
               endif
            endif
            dbdtlev = totplnkd(indlev+1,iband) - totplnkd(indlev,iband)
            dbdtlay = totplnkd(indlay+1,iband) - totplnkd(indlay,iband)
            planklayd(iplon, lay,iband) = totplnkd(indlay,iband) + tlayfrac * dbdtlay
            planklevd(iplon, lay,iband) = totplnkd(indlev,iband) + tlevfrac * dbdtlev
         enddo

!  For band 16, if radiative transfer will be performed on just
!  this band, use integrated Planck values up to 3250 cm-1.  
!  If radiative transfer will be performed across all 16 bands,
!  then include in the integrated Planck values for this band
!  contributions from 2600 cm-1 to infinity.
         iband = 16
         if (istart .eq. 16) then
            if (lay.eq.1) then
               dbdtlev = totplk16d( indbound+1) - totplk16d( indbound)
               plankbndd(iplon, iband) = semissd(iplon, iband) * &
                    (totplk16d( indbound) + tbndfrac * dbdtlev)
               if (idrv .eq. 1) then
                  dbdtlev = totplk16derivd( indbound+1) - totplk16derivd( indbound)
                  dplankbnd_dtd(iplon, iband) = semissd(iplon, iband) * &
                       (totplk16derivd(indbound) + tbndfrac * dbdtlev)
               endif
               dbdtlev = totplnkd(indlev0+1,iband)-totplnkd(indlev0,iband)
               planklevd(iplon, 0,iband) = totplk16d( indlev0) + &
                    t0frac * dbdtlev
            endif
            dbdtlev = totplk16d( indlev+1) - totplk16d( indlev)
            dbdtlay = totplk16d( indlay+1) - totplk16d( indlay)
            planklayd(iplon, lay,iband) = totplk16d( indlay) + tlayfrac * dbdtlay
            planklevd(iplon, lay,iband) = totplk16d( indlev) + tlevfrac * dbdtlev
         else
            if (lay.eq.1) then
               dbdtlev = totplnkd(indbound+1,iband) - totplnkd(indbound,iband)
               plankbndd(iplon, iband) = semissd(iplon, iband) * &
                    (totplnkd(indbound,iband) + tbndfrac * dbdtlev)
               if (idrv .eq. 1) then 
                  dbdtlev = totplnkderivd( indbound+1,iband) - totplnkderivd( indbound,iband)
                  dplankbnd_dtd(iplon, iband) = semissd(iplon, iband) * &
                       (totplnkderivd( indbound,iband) + tbndfrac * dbdtlev)
               endif
               dbdtlev = totplnkd(indlev0+1,iband)-totplnkd(indlev0,iband)
               planklevd(iplon, 0,iband) = totplnkd(indlev0,iband) + t0frac * dbdtlev
            endif
            dbdtlev = totplnkd(indlev+1,iband) - totplnkd(indlev,iband)
            dbdtlay = totplnkd(indlay+1,iband) - totplnkd(indlay,iband)
            planklayd(iplon, lay,iband) = totplnkd(indlay,iband) + tlayfrac * dbdtlay
            planklevd(iplon, lay,iband) = totplnkd(indlev,iband) + tlevfrac * dbdtlev
         endif

!  Find the two reference pressures on either side of the
!  layer pressure.  Store them in JP and JP1.  Store in FP the
!  fraction of the difference (in ln(pressure)) between these
!  two values that the layer pressure lies.
!         plog = alog(pavel(lay))
         plog = alog(pavel(iplon, lay))
         jp(iplon, lay) = int(36.  - 5*(plog+0.04 ))
         if (jp(iplon, lay) .lt. 1) then
            jp(iplon, lay) = 1
         elseif (jp(iplon, lay) .gt. 58) then
            jp(iplon, lay) = 58
         endif
         jp1 = jp(iplon, lay) + 1
         fp = 5.  *(preflogd(jp(iplon, lay)) - plog)

!  Determine, for each reference pressure (JP and JP1), which
!  reference temperature (these are different for each  
!  reference pressure) is nearest the layer temperature but does
!  not exceed it.  Store these indices in JT and JT1, resp.
!  Store in FT (resp. FT1) the fraction of the way between JT
!  (JT1) and the next highest reference temperature that the 
!  layer temperature falls.
         jt(iplon, lay) = int(3.  + (taveld(iplon, lay)-trefd(jp(iplon, lay)))/15. )
         if (jt(iplon, lay) .lt. 1) then
            jt(iplon, lay) = 1
         elseif (jt(iplon, lay) .gt. 4) then
            jt(iplon, lay) = 4
         endif
         ft = ((taveld(iplon, lay)-trefd(jp(iplon, lay)))/15. ) - float(jt(iplon, lay)-3)
         jt1(iplon, lay) = int(3.  + (taveld(iplon, lay)-trefd( jp1))/15. )
         if (jt1(iplon, lay) .lt. 1) then
            jt1(iplon, lay) = 1
         elseif (jt1(iplon, lay) .gt. 4) then
            jt1(iplon, lay) = 4
         endif
         ft1 = ((taveld(iplon, lay)-trefd(jp1))/15. ) - float(jt1(iplon, lay)-3)
         water = wv/lcoldry
         scalefac = pavel(iplon, lay) * stpfac / taveld(iplon, lay)

!  If the pressure is less than ~100mb, perform a different
!  set of species interpolations.
         if (plog .le. 4.56 ) go to 5300
         laytrop(iplon) =  laytrop(iplon) + 1

         forfac(iplon, lay) = scalefac / (1.+water)
         factor = (332.0 -taveld(iplon, lay))/36.0 
         indfor(iplon, lay) = min(2, max(1, int(factor)))
         forfrac(iplon, lay) = factor - float(indfor(iplon, lay))

!  Set up factors needed to separately include the water vapor
!  self-continuum in the calculation of absorption coefficient.
         selffac(iplon, lay) = water * forfac(iplon, lay)
         factor = (taveld(iplon, lay)-188.0 )/7.2 
         indself(iplon, lay) = min(9, max(1, int(factor)-7))
         selffrac(iplon, lay) = factor - float(indself(iplon, lay) + 7)

!  Set up factors needed to separately include the minor gases
!  in the calculation of absorption coefficient
         scaleminor(iplon, lay) = pavel(iplon, lay)/taveld(iplon, lay)
         scaleminorn2(iplon, lay) = (pavel(iplon, lay)/taveld(iplon, lay)) &
             *(wbroadd(iplon, lay)/(lcoldry+wv))
         factor = (taveld(iplon, lay)-180.8 )/7.2 
         indminor(iplon, lay) = min(18, max(1, int(factor)))
         minorfrac(iplon, lay) = factor - float(indminor(iplon, lay))

!  Setup reference ratio to be used in calculation of binary
!  species parameter in lower atmosphere.
         rat_h2oco2(iplon, lay)=chi_mlsd( 1,jp(iplon, lay))/chi_mlsd( 2,jp(iplon, lay))
         rat_h2oco2_1(iplon, lay)=chi_mlsd( 1,jp(iplon, lay)+1)/chi_mlsd( 2,jp(iplon, lay)+1)

         rat_h2oo3(iplon, lay)=chi_mlsd( 1,jp(iplon, lay))/chi_mlsd( 3,jp(iplon, lay))
         rat_h2oo3_1(iplon, lay)=chi_mlsd( 1,jp(iplon, lay)+1)/chi_mlsd( 3,jp(iplon, lay)+1)

         rat_h2on2o(iplon, lay)=chi_mlsd( 1,jp(iplon, lay))/chi_mlsd( 4,jp(iplon, lay))
         rat_h2on2o_1(iplon, lay)=chi_mlsd( 1,jp(iplon, lay)+1)/chi_mlsd( 4,jp(iplon, lay)+1)

         rat_h2och4(iplon, lay)=chi_mlsd( 1,jp(iplon, lay))/chi_mlsd( 6,jp(iplon, lay))
         rat_h2och4_1(iplon, lay)=chi_mlsd( 1,jp(iplon, lay)+1)/chi_mlsd( 6,jp(iplon, lay)+1)

         rat_n2oco2(iplon, lay)=chi_mlsd( 4,jp(iplon, lay))/chi_mlsd( 2,jp(iplon, lay))
         rat_n2oco2_1(iplon, lay)=chi_mlsd( 4,jp(iplon, lay)+1)/chi_mlsd( 2,jp(iplon, lay)+1)

!  Calculate needed column amounts.
         colh2o(iplon, lay) = 1.e-20  * colh2o(iplon, lay) * lcoldry
         colco2(iplon, lay) = 1.e-20  *  colco2(iplon, lay) * lcoldry
         colo3(iplon, lay) = 1.e-20  * colo3(iplon, lay) * lcoldry
         coln2o(iplon, lay) = 1.e-20  * coln2o(iplon, lay) * lcoldry
         colco(iplon, lay) = 1.e-20  * colco(iplon, lay) * lcoldry
         colch4(iplon, lay) = 1.e-20  * colch4(iplon, lay) * lcoldry
         colo2(iplon, lay) = 1.e-20  * colo2(iplon, lay) * lcoldry
         if (colco2(iplon, lay) .eq. 0. ) colco2(iplon, lay) = 1.e-32  * lcoldry
         if (colo3(iplon, lay) .eq. 0. ) colo3(iplon, lay) = 1.e-32  * lcoldry
         if (coln2o(iplon, lay) .eq. 0. ) coln2o(iplon, lay) = 1.e-32  * lcoldry
         if (colco(iplon, lay) .eq. 0. ) colco(iplon, lay) = 1.e-32  * lcoldry
         if (colch4(iplon, lay) .eq. 0. ) colch4(iplon, lay) = 1.e-32  * lcoldry
         colbrd(iplon, lay) = 1.e-20  * wbroadd(iplon, lay)
         go to 5400

!  Above laytrop.
 5300    continue

         forfac(iplon, lay) = scalefac / (1.+water)
         factor = (taveld(iplon, lay)-188.0 )/36.0 
         indfor(iplon, lay) = 3
         forfrac(iplon, lay) = factor - 1.0 

!  Set up factors needed to separately include the water vapor
!  self-continuum in the calculation of absorption coefficient.
         selffac(iplon, lay) = water * forfac(iplon, lay)

!  Set up factors needed to separately include the minor gases
!  in the calculation of absorption coefficient
         scaleminor(iplon, lay) = pavel(iplon, lay)/taveld(iplon, lay)         
         scaleminorn2(iplon, lay) = (pavel(iplon, lay)/taveld(iplon, lay)) &
             * (wbroadd(iplon, lay)/(coldry(iplon, lay)+wv))
         factor = (taveld(iplon, lay)-180.8 )/7.2 
         indminor(iplon, lay) = min(18, max(1, int(factor)))
         minorfrac(iplon, lay) = factor - float(indminor(iplon, lay))

!  Setup reference ratio to be used in calculation of binary
!  species parameter in upper atmosphere.
         rat_h2oco2(iplon, lay)=chi_mlsd( 1,jp(iplon, lay))/chi_mlsd( 2,jp(iplon, lay))
         rat_h2oco2_1(iplon, lay)=chi_mlsd( 1,jp(iplon, lay)+1)/chi_mlsd( 2,jp(iplon, lay)+1)         

         rat_o3co2(iplon, lay)=chi_mlsd( 3,jp(iplon, lay))/chi_mlsd( 2,jp(iplon, lay))
         rat_o3co2_1(iplon, lay)=chi_mlsd( 3,jp(iplon, lay)+1)/chi_mlsd( 2,jp(iplon, lay)+1)         

!  Calculate needed column amounts.
         colh2o(iplon, lay) = 1.e-20  * colh2o(iplon, lay) * lcoldry
         colco2(iplon, lay) = 1.e-20  *  colco2(iplon, lay) * lcoldry
         colo3(iplon, lay) = 1.e-20  * colo3(iplon, lay) * lcoldry
         coln2o(iplon, lay) = 1.e-20  * coln2o(iplon, lay) * lcoldry
         colco(iplon, lay) = 1.e-20  * colco(iplon, lay) * lcoldry
         colch4(iplon, lay) = 1.e-20  * colch4(iplon, lay) * lcoldry
         colo2(iplon, lay) = 1.e-20  * colo2(iplon, lay) * lcoldry
         if (colco2(iplon, lay) .eq. 0. ) colco2(iplon, lay) = 1.e-32  * lcoldry
         if (colo3(iplon, lay) .eq. 0. ) colo3(iplon, lay) = 1.e-32  * lcoldry
         if (coln2o(iplon, lay) .eq. 0. ) coln2o(iplon, lay) = 1.e-32  * lcoldry
         if (colco(iplon, lay)  .eq. 0. ) colco(iplon, lay) = 1.e-32  * lcoldry
         if (colch4(iplon, lay) .eq. 0. ) colch4(iplon, lay) = 1.e-32  * lcoldry
         colbrd(iplon, lay) = 1.e-20  * wbroadd(iplon, lay)
 5400    continue

!  We have now isolated the layer ln pressure and temperature,
!  between two reference pressures and two reference temperatures 
!  (for each reference pressure).  We multiply the pressure 
!  fraction FP with the appropriate temperature fractions to get 
!  the factors that will be needed for the interpolation that yields
!  the optical depths (performed in routines TAUGBn for band n).`

         compfp = 1. - fp
         fac10(iplon, lay) = compfp * ft
         fac00(iplon, lay) = compfp * (1.  - ft)
         fac11(iplon, lay) = fp * ft1
         fac01(iplon, lay) = fp * (1.  - ft1)

!  Rescale selffac and forfac for use in taumol
         selffac(iplon, lay) = colh2o(iplon, lay)*selffac(iplon, lay)
         forfac(iplon, lay) = colh2o(iplon, lay)*forfac(iplon, lay)

	
! End layer loop
      enddo
     
#ifdef _CUDA
    endif
#else
    end do
#endif
      end subroutine setcoefg

end module gpu_rrtmg_lw_setcoef