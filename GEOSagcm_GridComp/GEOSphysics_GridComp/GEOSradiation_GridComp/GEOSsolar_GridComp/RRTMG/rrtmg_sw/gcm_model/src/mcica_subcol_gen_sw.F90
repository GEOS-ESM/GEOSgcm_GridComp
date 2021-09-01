module mcica_subcol_gen_sw

   ! Purpose: Create stochastic arrays for cloud physical properties.
   ! Input gridcolumn cloud profiles: cloud fraction and in-cloud ice and liquid
   ! water paths. Output will be stochastic subcolumn arrays of these variables.

   ! PMN 2021/08: This version makes the binary (cldf_stoch in {0,1}) nature of
   ! cloud output explicit by changing to a logical cldy_stoch in {true,false},
   ! which may permit some exterior efficiency improvements.

   use cloud_condensate_inhomogeneity, only: &
      condensate_inhomogeneous, zcw_lookup

   implicit none
   private

   public :: mcica_sw

contains

   !---------------------------------------------------------------------------------------
   subroutine mcica_sw( &
      dncol, ncol, nsubcol, nlay, &
      zmid, alat, doy, &
      play, cldfrac, ciwp, clwp, cwp_tiny, &
      cldy_stoch, ciwp_stoch, clwp_stoch)
   !---------------------------------------------------------------------------------------
   !
   ! Original code: Based on Raisanen et al., QJRMS, 2004.
   ! Original Contact: Cecile Hannay (hannay@ucar.edu)
   !
   ! Modifications:
   ! ...
   ! Peter Norris, GMAO, Apr-Aug 2021:
   !   Tidy up code and comments
   !   Keep only exponential (generalized) overlap.
   !   Abstract and reorder indices for CPU efficiency.
   !
   !   Given a profile of cloud fraction, cloud water and cloud ice, produce a set of
   ! subcolumns. Each subcolumn has cloud fraction in {0,1} at each level.
   !   If not condensate_inhomogeneous(), each layer within each subcolumn has uniform
   ! cloud ice and cloud liquid concentration. If condensate_inhomogeneous(), each layer
   ! has horizontal condensate variability.
   !   The ensemble of subcolumns statistically reproduces the cloud fraction within each
   ! layer and its prescribed vertical overlap, and, if condensate_inhomogeneous(), the
   ! PDF of cloud liquid and ice within each layer and its prescribed vertical correlation.
   ! Whether there is condensate inhomogeneity, and the condensate PDF if so, are set in
   ! initialize_inhomogeneity().
   !
   ! Overlap assumption:
   ! Exponential (generalized) overlap (Raisanen et al., 2004, Pincus et al., 2005) using
   ! correlations alpha and rcorr based on decorrelation length scales from Oreopoulos et
   ! al., 2012.
   !
   !---------------------------------------------------------------------------------------

      integer, intent(in) :: dncol                 ! Dimensioned number of gridcols
      integer, intent(in) :: ncol                  ! Actual number of gridcols
      integer, intent(in) :: nsubcol               ! #Subcols to generate / gridcol
      integer, intent(in) :: nlay                  ! Number of model layers

      real,    intent(in) :: zmid    (nlay,dncol)  ! Height of midpoints [m]
      real,    intent(in) :: alat         (dncol)  ! Latitude of gridcolumn [radians]
      integer, intent(in) :: doy                   ! Day of year
      real,    intent(in) :: play    (nlay,dncol)  ! Layer pressures [hPa]
      real,    intent(in) :: cldfrac (nlay,dncol)  ! Layer cloud fraction [0., 1.]

      ! The units of these in-cloud water paths are not specified, but they should
      ! be the same for both liquid water and ice. cwp_tiny below is assumed to be
      ! in those isame units and the output stochastic water paths as well.

      real,    intent(in) :: ciwp    (nlay,dncol)  ! In-cloud ice water path
      real,    intent(in) :: clwp    (nlay,dncol)  ! In-cloud liquid water path

      ! Tiny threshold value for cloud water path, at or below which each of
      ! ice and liquid water parths are separately reset to zero. If both are
      ! thus reset, then the subcolumn gridcell is reset to .not.cloudy. The
      ! idea here is mainly efficiency, since external processing of cloudy
      ! cells usually takes longer. Units are assumed same as ciwp and clwp.

      real,    intent(in) :: cwp_tiny

      ! output subcolumns
      ! (units of water paths are the same as for inputs ciwp and clwp)
      logical, intent(out) :: cldy_stoch (nlay,nsubcol,dncol)  ! Cloudy or not?
      real,    intent(out) :: ciwp_stoch (nlay,nsubcol,dncol)  ! In-cloud ice water path
      real,    intent(out) :: clwp_stoch (nlay,nsubcol,dncol)  ! In-cloud liq water path

      ! ----- Locals -----

      ! decorrelation length scales for cldfrac and condensate
      real, dimension(dncol) :: adl, rdl

      ! inter-layer correlations for cldfrac and condensate
      real, dimension(nlay,dncol) :: alpha, rcorr

      ! seeds and random number (rng_kiss)
      integer :: seed1, seed2, seed3, seed4
      real :: rand_num

      ! random number arrays used for overlap
      real, dimension(nlay,nsubcol,dncol) :: &
         cdf1, &  ! for cloud presence
         cdf2, &  ! auxilliary
         cdf3     ! for cloud condensate

      ! other locals
      real :: cfs, sigma_qcw, zcw
      logical :: cond_inhomo, ciwp_negligible, clwp_negligible

      ! indices
      integer :: icol, isubcol, ilay

      ! save for speed
      cond_inhomo = condensate_inhomogeneous()

      ! -----------------------------------
      ! compute decorrelation length scales
      ! -----------------------------------

      ! for cloud presence
      call correlation_length(dncol, ncol, &
         1.4315, 2.1219, -25.584, 7., doy, alat, adl)

      ! for condensate
      if (cond_inhomo) then
         call correlation_length(dncol, ncol, &
            0.72192, 0.78996, 40.404, 8.5, doy, alat, rdl)
      endif

      ! --------------------------
      ! outer loop over gridcolumn
      ! --------------------------
      do icol = 1,ncol

         ! ------------------------------------
         ! exponential inter-layer correlations
         ! ------------------------------------
         do ilay = 2,nlay
            alpha(ilay,icol) = exp( -(zmid(ilay,icol)-zmid(ilay-1,icol)) / adl(icol) )
         enddo
         if (cond_inhomo) then
            do ilay = 2,nlay
               rcorr(ilay,icol) = exp( -(zmid(ilay,icol)-zmid(ilay-1,icol)) / rdl(icol) )
            enddo
         endif

         do isubcol = 1,nsubcol

            seed1 = (play(1,icol)*100. - int(play(1,icol)*100.)) * 1000000000 + isubcol * 11
            seed3 = (play(3,icol)*100. - int(play(3,icol)*100.)) * 1000000000 + isubcol * 13
            seed2 = seed1 + isubcol
            seed4 = seed3 - isubcol

            do ilay = 1,nlay
               call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
               cdf1(ilay,isubcol,icol) = rand_num
               call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
               cdf2(ilay,isubcol,icol) = rand_num
            enddo

            do ilay = 2,nlay
               if (cdf2(ilay,isubcol,icol) < alpha(ilay,icol)) then
                  cdf1(ilay,isubcol,icol) = cdf1(ilay-1,isubcol,icol) 
               endif
            enddo
         enddo

      enddo
    
      if (cond_inhomo) then
        
         do icol = 1,ncol

            do isubcol = 1,nsubcol
               seed1 = (play(1,icol)*100. - int(play(1,icol)*100.)) * 1000000000 + isubcol * 11
               seed3 = (play(3,icol)*100. - int(play(3,icol)*100.)) * 1000000000 + isubcol * 13
               seed2 = seed1 + isubcol
               seed4 = seed3 - isubcol
               do ilay = 1,nlay
                  call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
                  cdf2(ilay,isubcol,icol) = rand_num
                  call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
                  cdf3(ilay,isubcol,icol) = rand_num
               enddo

               do ilay = 2,nlay
                  if (cdf2(ilay,isubcol,icol) < rcorr(ilay,icol)) then
                     cdf3(ilay,isubcol,icol) = cdf3(ilay-1,isubcol,icol)
                  endif
               enddo
            enddo

         enddo

      endif

      ! -------------------
      ! generate subcolumns
      ! -------------------

      do icol = 1,ncol
         do isubcol = 1,nsubcol
            do ilay = 1,nlay
               cfs = cldfrac(ilay,icol)

               if (cdf1(ilay,isubcol,icol) >= 1. - cfs) then

                  ! a cloudy subcolumn/layer with inhomogeneous condensate assignment
                  ! note: cdf1 from rand_num in (0,1) so:
                  !       cldfrac == 0. can never get here;
                  !       cldfrac == 1. always comes here.

                  cldy_stoch(ilay,isubcol,icol) = .true.

                  if (cond_inhomo) then

                     ! Cloud fraction sets level of inhomogeneity
                     if (cfs > 0.99) then
                        sigma_qcw = 0.5
                     elseif (cfs > 0.9) then
                        sigma_qcw = 0.71
                     else
                        sigma_qcw = 1.0
                     endif

                     ! horizontally variable clouds
                     zcw = zcw_lookup(cdf3(ilay,isubcol,icol),sigma_qcw)
                     ciwp_stoch(ilay,isubcol,icol) = ciwp(ilay,icol) * zcw
                     clwp_stoch(ilay,isubcol,icol) = clwp(ilay,icol) * zcw

                  else

                     ! horizontally homogeneous clouds
                     ciwp_stoch(ilay,isubcol,icol) = ciwp(ilay,icol)
                     clwp_stoch(ilay,isubcol,icol) = clwp(ilay,icol)

                  endif

                  ! reset negligible water paths to zero (see comment in introduction)
                  ciwp_negligible = (ciwp_stoch(ilay,isubcol,icol) <= cwp_tiny)
                  if (ciwp_negligible) ciwp_stoch(ilay,isubcol,icol) = 0.
                  clwp_negligible = (clwp_stoch(ilay,isubcol,icol) <= cwp_tiny)
                  if (clwp_negligible) clwp_stoch(ilay,isubcol,icol) = 0.

                  ! reset cloudy status if both negligible
                  if (ciwp_negligible .and. clwp_negligible) &
                     cldy_stoch(ilay,isubcol,icol) = .false.

               else

                  ! a clear subcolumn/layer
                  ! note: per comments above:
                  !       cldfrac == 0. always comes here;
                  !       cldfrac == 1. never comes here.

                  ! a clear subcolumn
                  cldy_stoch(ilay,isubcol,icol) = .false.
                  ciwp_stoch(ilay,isubcol,icol) = 0.
                  clwp_stoch(ilay,isubcol,icol) = 0.

               endif
            enddo
         enddo
      enddo

   end subroutine mcica_sw

   !-------------------------------------------
   subroutine correlation_length(dncol, ncol, &
      am1, am2, am4, amr, doy, alat, clength)
   !-------------------------------------------
      integer, intent(in)  :: dncol               ! Dimensioned number of gridcols
      integer, intent(in)  :: ncol                ! Actual number of gridcols
      real,    intent(in)  :: am1, am2, am4, amr  ! input parameters
      integer, intent(in)  :: doy                 ! Day of year
      real,    intent(in)  :: alat    (dncol)     ! Latitude of gridcolumn [radians]
      real,    intent(out) :: clength (dncol)     ! Correlation length [m]

      real, parameter :: r2d = 180.d0 / 3.14159265358979323846d0

      integer :: icol
      real :: am3

      if (doy > 181) then
         am3 = -4.*amr/365.*(doy-272)
      else
         am3 =  4.*amr/365.*(doy- 91)
      endif

      do icol = 1,ncol
         clength(icol) = (am1+am2*exp(-(alat(icol)*r2d-am3)**2/am4**2))*1.e3
      enddo

   end subroutine correlation_length


   !-------------------------------------------------------
   subroutine rng_kiss(seed1, seed2, seed3, seed4, ran_num)
   !-------------------------------------------------------
   ! See note below: get ran_num in (0,1)

      real, intent(out) :: ran_num
      integer, intent(inout) :: seed1, seed2, seed3, seed4
      integer :: m, k, n, kiss

      ! inline function
      m(k,n) = ieor(k,ishft(k,n))

      seed1 = 69069 * seed1 + 1327217885
      seed2 = m(m(m(seed2,13),-17),5)
      seed3 = 18000 * iand(seed3,65535) + ishft(seed3,-16)
      seed4 = 30903 * iand(seed4,65535) + ishft(seed4,-16)
      kiss = seed1 + seed2 + ishft(seed3,16) + seed4
      ran_num = kiss * 2.328306e-10 + 0.5

   end subroutine rng_kiss

end module mcica_subcol_gen_sw
