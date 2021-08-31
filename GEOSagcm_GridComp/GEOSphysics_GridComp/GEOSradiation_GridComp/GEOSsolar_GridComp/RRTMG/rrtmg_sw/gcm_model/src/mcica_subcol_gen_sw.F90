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

   ! tiny threshold value for cloud water path, at or below which each of
   ! ice and liquid water parths are separately reset to zero. If both are
   ! thus reset, then the subcolumn gridcell is reset to .not.cloudy. The
   ! idea here is mainly efficiency, since external processing of cloudy
   ! cells usually takes longer.
   real, parameter :: cwp_tiny = 1.e-20
!pmn? units

   real, parameter :: r2d = 180.d0 / 3.14159265358979323846d0

   public :: mcica_sw
   public :: cwp_tiny

contains

   !---------------------------------------------------------------------------------------
   subroutine mcica_sw( &
      pncol, ncol, nsubcol, nlay, &
      zmid, alat, doy, &
      play, cldfrac, ciwp, clwp, &
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

      integer, intent(in) :: pncol                 ! Dimensioned number of gridcols
      integer, intent(in) :: ncol                  ! Actual number of gridcols
      integer, intent(in) :: nsubcol               ! #Subcols to generate / gridcol
      integer, intent(in) :: nlay                  ! Number of model layers
      real,    intent(in) :: zmid    (nlay,pncol)  ! Hgt of midpoints [m]
      real,    intent(in) :: alat         (pncol)  ! Latitude of gridcolumn [radians]
      integer, intent(in) :: doy                   ! Day of year
      real,    intent(in) :: play    (nlay,pncol)  ! Layer pressures [hPa]
!pmn: original comment had [Pa] but [hPa] is being passed in !!!
      real,    intent(in) :: cldfrac (nlay,pncol)  ! Layer cloud fraction 
      real,    intent(in) :: ciwp    (nlay,pncol)  ! In-cloud ice water path (g/m2)?
      real,    intent(in) :: clwp    (nlay,pncol)  ! In-cloud liquid water path (g/m2)?

      ! output subcolumns
      ! (units of water paths are the same as for inputs ciwp and clwp)
      logical, intent(out) :: cldy_stoch (nlay,nsubcol,pncol)  ! Cloudy or not?
      real,    intent(out) :: ciwp_stoch (nlay,nsubcol,pncol)  ! In-cloud ice water path
      real,    intent(out) :: clwp_stoch (nlay,nsubcol,pncol)  ! In-cloud liq water path

      ! ----- Locals -----

      ! decorrelation length scales for cldfrac and condensate
      real, dimension(pncol) :: adl, rdl

      ! inter-layer correlations for cldfrac and condensate
      real, dimension(nlay,pncol) :: alpha, rcorr

      ! related to decorrelation lengths
      real :: am1, am2, am3, am4, amr

      ! related to random number and seed 
      integer :: seed1, seed2, seed3, seed4  ! seed (rng_kiss)
      real :: rand_num                       ! random number (rng_kiss)

      ! random number arrays used for overlap
      real, dimension(nlay,nsubcol,pncol) :: &
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

      cdf1 = 0.
      cdf2 = 0.
      cdf3 = 0.
      alpha = 0.
      rcorr = 0.

      ! -----------------------------------
      ! compute decorrelation length scales
      ! -----------------------------------

      ! cloud presence decorrelation length scale
      am1 = 1.4315
      am2 = 2.1219
      am4 = -25.584
      amr = 7.
      if (doy > 181) then
         am3 = -4.*amr/365*(doy-272)
      else
         am3 = 4.*amr/365.*(doy-91)
      endif
      do icol = 1,ncol
         adl(icol) = (am1+am2*exp(-(alat(icol)*r2d-am3)**2/am4**2))*1.e3  ! [m]
      end do

      ! condensate decorrelation length scale
      am1 = 0.72192
      am2 = 0.78996
      am4 = 40.404
      amr = 8.5
      if (doy > 181) then
         am3 = -4.*amr/365*(doy-272)
      else
         am3 = 4.*amr/365.*(doy-91)
      endif
      do icol = 1,ncol
         rdl(icol) = (am1+am2*exp(-(alat(icol)*r2d-am3)**2/am4**2))*1.e3  ! [m]
      end do
   
      do icol = 1,ncol
         alpha(1,icol) = 0.
         do ilay = 2,nlay
            alpha(ilay,icol) = exp(-(zmid(ilay,icol) - zmid(ilay-1,icol)) / adl(icol))
         end do
      end do
     
! pmn eventually reorder but not now cos will break zero-diff
! preserve zero-diff as long as possible so work on other parts of code first

      do icol = 1,ncol
         do isubcol = 1,nsubcol

!        do ilay = 1,nlay
! play here is actually in hPa, so ~ 950-1050 at sfc.
! difference clearly [0,1) (a fractional hPa).
! so after mult [0,100000000) cf range of i4: [-2147483648,2147483647],
! which is well within range --- have 10 digits use only 9
!           seed1 = (play(1,icol) - int(play(1,icol))) * 100000000 - ilay
!           seed2 = (play(2,icol) - int(play(2,icol))) * 100000000 + ilay
!           seed3 = (play(3,icol) - int(play(3,icol))) * 100000000 + ilay * 6.2
!           seed4 = (play(4,icol) - int(play(4,icol))) * 100000000           
!pmn: 8 zeros here cf 9 for LW!!!

            seed1 = (play(1,icol)*100. - int(play(1,icol)*100.)) * 1000000000 + isubcol * 11
            seed3 = (play(3,icol)*100. - int(play(3,icol)*100.)) * 1000000000 + isubcol * 13
            seed2 = seed1 + isubcol
            seed4 = seed3 - isubcol

!           do isubcol = 1,nsubcol
            do ilay = 1,nlay
               call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
               cdf1(ilay,isubcol,icol) = rand_num 
               call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
               cdf2(ilay,isubcol,icol) = rand_num 
            end do
         end do
      end do
    
      do icol = 1,ncol
         do isubcol = 1,nsubcol
            do ilay = 2,nlay
               if (cdf2(ilay,isubcol,icol) < alpha(ilay,icol)) then
                  cdf1(ilay,isubcol,icol) = cdf1(ilay-1,isubcol,icol) 
               end if
            end do
         end do
      end do
    
      if (cond_inhomo) then
        
         do icol = 1,ncol
            rcorr(1,icol) = 0.
            do ilay = 2,nlay
               rcorr(ilay,icol) = exp(-(zmid(ilay,icol) - zmid(ilay-1,icol)) / rdl(icol))
            end do
         end do
        
         do icol = 1,ncol
            do isubcol = 1,nsubcol
!           do ilay = 1,nlay
!              seed1 = (play(1,icol) - int(play(1,icol))) * 100000000 - ilay
!              seed2 = (play(2,icol) - int(play(2,icol))) * 100000000 + ilay
!              seed3 = (play(3,icol) - int(play(3,icol))) * 100000000 + ilay * 6.2
!              seed4 = (play(4,icol) - int(play(4,icol))) * 100000000           
               seed1 = (play(1,icol)*100. - int(play(1,icol)*100.)) * 1000000000 + isubcol * 11
               seed3 = (play(3,icol)*100. - int(play(3,icol)*100.)) * 1000000000 + isubcol * 13
               seed2 = seed1 + isubcol
               seed4 = seed3 - isubcol
!              do isubcol = 1,nsubcol
               do ilay = 1,nlay
                  call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
                  cdf2(ilay,isubcol,icol) = rand_num 
                  call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
                  cdf3(ilay,isubcol,icol) = rand_num 
               end do
            end do
         end do

         do icol = 1,ncol
            do isubcol = 1,nsubcol
               do ilay = 2,nlay
                  if (cdf2(ilay,isubcol,icol) < rcorr(ilay,icol)) then
                     cdf3(ilay,isubcol,icol) = cdf3(ilay-1,isubcol,icol)
                  end if
               end do
            end do
         end do
      end if

      ! generate clear or cloudy cells of subcolumns

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
                     clwp_stoch(ilay,isubcol,icol) = clwp(ilay,icol) * zcw
                     ciwp_stoch(ilay,isubcol,icol) = ciwp(ilay,icol) * zcw

                  else

                     ! horizontally homogeneous clouds
                     clwp_stoch(ilay,isubcol,icol) = clwp(ilay,icol)
                     ciwp_stoch(ilay,isubcol,icol) = ciwp(ilay,icol)

                  end if
                
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
                  clwp_stoch(ilay,isubcol,icol) = 0. 
                  ciwp_stoch(ilay,isubcol,icol) = 0. 

               endif
            enddo
         enddo
      enddo

   end subroutine mcica_sw

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
