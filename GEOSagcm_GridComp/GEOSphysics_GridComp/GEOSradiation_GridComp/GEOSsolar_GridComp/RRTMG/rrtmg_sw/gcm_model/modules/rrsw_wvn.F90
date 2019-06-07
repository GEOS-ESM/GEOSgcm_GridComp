      module rrsw_wvn

      !use parkind, only : im => kind , rb => kind 
      use parrrsw, only : nbndsw, mg, ngptsw, jpb1, jpb2

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw spectral information

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ng     :  integer: Number of original g-intervals in each spectral band
! nspa   :  integer: 
! nspb   :  integer: 
!wavenum1:  real   : Spectral band lower boundary in wavenumbers
!wavenum2:  real   : Spectral band upper boundary in wavenumbers
! delwave:  real   : Spectral band width in wavenumbers
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
!                    (224 total) into reduced set of g-intervals 
!                    (112 total)
!------------------------------------------------------------------

      integer  :: ng(jpb1:jpb2)
      integer  :: nspa(jpb1:jpb2)
      integer  :: nspb(jpb1:jpb2)

      real  :: wavenum1(jpb1:jpb2)
      real  :: wavenum2(jpb1:jpb2)
      real  :: delwave(jpb1:jpb2)
      integer :: icxa(jpb1:jpb2)

      integer  :: ngc(nbndsw)
      integer  :: ngs(nbndsw)
      integer  :: ngn(ngptsw)
      integer  :: ngb(ngptsw)
      integer  :: ngm(nbndsw*mg)

      real  :: wt(mg)
      real  :: rwgt(nbndsw*mg)

      end module rrsw_wvn
