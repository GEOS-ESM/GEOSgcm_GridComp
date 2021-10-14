      module rrlw_wvn

       !use parkind, only : im => kind , rb => kind 
      use parrrtm, only : nbndlw, mg, ngptlw, maxinpx

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw spectral information

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ng     :  integer: Number of original g-intervals in each spectral band
! nspa   :  integer: For the lower atmosphere, the number of reference
!                    atmospheres that are stored for each spectral band
!                    per pressure level and temperature.  Each of these
!                    atmospheres has different relative amounts of the 
!                    key species for the band (i.e. different binary
!                    species parameters).
! nspb   :  integer: Same as nspa for the upper atmosphere
!wavenum1:  real   : Spectral band lower boundary in wavenumbers
!wavenum2:  real   : Spectral band upper boundary in wavenumbers
! delwave:  real   : Spectral band width in wavenumbers
! totplnk:  real   : Integrated Planck value for each band; (band 16
!                    includes total from 2600 cm-1 to infinity)
!                    Used for calculation across total spectrum
!totplk16:  real   : Integrated Planck value for band 16 (2600-3250 cm-1)
!                    Used for calculation in band 16 only if 
!                    individual band output requested
!totplnkderiv: real: Integrated Planck function derivative with respect
!                    to temperature for each band; (band 16
!                    includes total from 2600 cm-1 to infinity)
!                    Used for calculation across total spectrum
!totplk16deriv:real: Integrated Planck function derivative with respect
!                    to temperature for band 16 (2600-3250 cm-1)
!                    Used for calculation in band 16 only if 
!                    individual band output requested
!
! ngc    :  integer: The number of new g-intervals in each band
! ngs    :  integer: The cumulative sum of new g-intervals for each band
! ngm    :  integer: The index of each new g-interval relative to the
!                    original 16 g-intervals in each band
! ngn    :  integer: The number of original g-intervals that are 
!                    combined to make each new g-intervals in each band
! ngb    :  integer: The band index for each new g-interval
! wt     :  real   : RRTM weights for the original 16 g-intervals
! rwgt   :  real   : Weights for combining original 16 g-intervals 
!                    (256 total) into reduced set of g-intervals 
!                    (140 total)
! nxmol  :  integer: Number of cross-section molecules
! ixindx :  integer: Flag for active cross-sections in calculation
!------------------------------------------------------------------

      integer  :: ng(nbndlw)
      integer  :: nspa(nbndlw)
      integer  :: nspb(nbndlw)

!     real  :: wavenum1(nbndlw)
!     real  :: wavenum2(nbndlw)
!     real  :: delwave(nbndlw)


      ! lower band limit [cm-1]
      real, parameter :: wavenum1 (nbndlw) = &
         [  10.,  350.,  500.,  630.,  700.,  820.,  980., 1080., &
          1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600.  ]
      ! upper band limit [cm-1]
      real, parameter :: wavenum2 (nbndlw) = &
         [ 350.,  500.,  630.,  700.,  820.,  980., 1080., 1180., &
          1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250.  ]
      ! band width [cm-1]
      real, parameter :: delwave  (nbndlw) = wavenum2 - wavenum1


      real  :: totplnk(181,nbndlw)
      real  :: totplk16(181)

      real  :: totplnkderiv(181,nbndlw)
      real  :: totplk16deriv(181)

      integer  :: ngc(nbndlw)
      integer  :: ngs(nbndlw)
      integer  :: ngn(ngptlw)
      integer  :: ngb(ngptlw)
      integer  :: ngm(nbndlw*mg)

      real  :: wt(mg)
      real  :: rwgt(nbndlw*mg)

      integer  :: nxmol
      integer  :: ixindx(maxinpx)

      end module rrlw_wvn
