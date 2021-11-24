!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

module cloud_subcol_gen

   ! Purpose: Create stochastic arrays for cloud physical properties.
   ! Input gridcolumn cloud profiles: cloud fraction and in-cloud ice and liquid
   ! water paths. Output will be stochastic subcolumn arrays of these variables.

   ! PMN 2021/06: This version makes the binary (cldf_stoch in {0,1}) nature of
   ! cloud output explicit by changing to a logical cldy_stoch in {true,false},
   ! which may permit some exterior efficiency improvements.

   use cloud_condensate_inhomogeneity, only: &
      condensate_inhomogeneous, zcw_lookup
   use iso_fortran_env, only : error_unit

   implicit none

   ! state is saved outside of scope
   save

   ! all fields and methods hidden unless declared public
   private

   ! Original correlation length parameters. 
   ! These are the full precision versions of those presented in
   ! Oreopoulos et al., Atmos. Chem. Phys. (2012) and SHOULD NOT BE
   ! CHANGED. They are for historical reference only. If you want
   ! to change the default values, change the def_ parameters below.
   ! Cloud presence:
   real, parameter :: Oreo12_aam1  = 1.4315
   real, parameter :: Oreo12_aam2  = 2.1219
   real, parameter :: Oreo12_aam30 = 7.
   real, parameter :: Oreo12_aam4  = -25.584
   ! Cloud condensate:
   real, parameter :: Oreo12_ram1  = 0.72192
   real, parameter :: Oreo12_ram2  = 0.78996
   real, parameter :: Oreo12_ram30 = 8.5
   real, parameter :: Oreo12_ram4  = 40.404

   ! Default correlation length parameters. 
   ! These are public so they can be accessed externally as defaults
   ! in MAPL_GetResource(). These can be changed from the historical
   ! Oreo12 values above if better defaults are found through tuning.
   ! Cloud presence:
   real, public, parameter :: def_aam1  = Oreo12_aam1
   real, public, parameter :: def_aam2  = Oreo12_aam2
   real, public, parameter :: def_aam30 = Oreo12_aam30
   real, public, parameter :: def_aam4  = Oreo12_aam4
   ! Cloud condensate:
   real, public, parameter :: def_ram1  = Oreo12_ram1
   real, public, parameter :: def_ram2  = Oreo12_ram2
   real, public, parameter :: def_ram30 = Oreo12_ram30
   real, public, parameter :: def_ram4  = Oreo12_ram4

   ! Actual correlation length parameters used.
   ! These are set to the default values above initially.
   ! To use different values, from Radiation GC Initialize do
   !     use cloud_subcol_gen, only : &
   !       initialize_cloud_subcol_gen, def_aam1, ...
   !     real :: aam1, ...
   !     call MAPL_GetResource(MAPL,aam1,LABEL="ADL_AM1:",default=def_aam1,__RC__)
   !     ...
   !     call initialize_cloud_subcol_gen(adl_am1=aam1, ...)
   ! This is only necessary to CHANGE the defaults above. 
   ! Cloud presence:
   real :: aam1  = def_aam1
   real :: aam2  = def_aam2
   real :: aam30 = def_aam30
   real :: aam4  = def_aam4
   ! Cloud condensate:
   real :: ram1  = def_ram1
   real :: ram2  = def_ram2
   real :: ram30 = def_ram30
   real :: ram4  = def_ram4

   ! public interface
   public :: initialize_cloud_subcol_gen  ! set non-default correlation lengths
   public :: generate_stochastic_clouds   ! generate subcolumns
   public :: clearCounts_threeBand        ! L|M|H|T cloud fractions

contains

   !---------------------------------------------------------------------------------------
   subroutine initialize_cloud_subcol_gen( &
      adl_am1, adl_am2, adl_am30, adl_am4, &
      rdl_am1, rdl_am2, rdl_am30, rdl_am4)
   !---------------------------------------------------------------------------------------
   ! Use if want to set non-default correlation lengths.
   !---------------------------------------------------------------------------------------

      ! correlation length parameters
      real, intent(in), optional :: adl_am1, adl_am2, adl_am30, adl_am4  ! cloud presence
      real, intent(in), optional :: rdl_am1, rdl_am2, rdl_am30, rdl_am4  ! cloud condensate

      ! optionally reset correlation length parameters from module level defaults
      if (present(adl_am1 )) aam1  = adl_am1
      if (present(adl_am2 )) aam2  = adl_am2
      if (present(adl_am30)) aam30 = adl_am30
      if (present(adl_am4 )) aam4  = adl_am4
      if (present(rdl_am1 )) ram1  = rdl_am1
      if (present(rdl_am2 )) ram2  = rdl_am2
      if (present(rdl_am30)) ram30 = rdl_am30
      if (present(rdl_am4 )) ram4  = rdl_am4

   end subroutine initialize_cloud_subcol_gen

   !---------------------------------------------------------------------------------------
   subroutine generate_stochastic_clouds( &
      dncol, ncol, nsubcol, nlay, &
      zmid, alat, doy, &
      play, cldfrac, ciwp, clwp, cwp_tiny, &
      cldy_stoch, ciwp_stoch, clwp_stoch, &
      seed_order)
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
   ! Peter Norris, GMAO, Nov 2021:
   !   Added optional seed_order (sse below).
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
   ! Optional seed order:
   ! This is an array of 4 integers which MUST be a permutation of [1,2,3,4]. If absent,
   ! [1,2,3,4] is assumed. It is used in the code to permute the pseed(4) array as part
   ! of the initial seeding on each gridcolumn's PRNG stream. As such, it provides a
   ! way of generating a different set of random numbers for say LW and SW at the same
   ! timetep, when the layer pressures <play> used for seed generation are the same.
   ! For example, LW can be run with [1,2,3,4] and SW with [4,2,1,3]. 
   ! NOTE: Using the same seed_order for LW and SW still does not mean using the "same
   ! subcolumn cloud ensemble" because the number of subcolumns (g-points) in LW and SW
   ! differ. Also, of course, McICA .ne. ICA, meaning that each subcolumn in McICA gets
   ! a different g-point, and those g-points have different meanings in LW and SW. So
   ! the same_seed order in LW and SW would only use an identical cloud field in the
   ! case of full ICA (each subcolumn sees every g-point) and with the same nsubcol,
   ! neither of which is true. NEVERTHELESS, if the user is worried that using the same 
   ! seeds for LW and SW will create an unwanted correlation of some sort, they can use
   ! a different seed_order for each. 
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
      ! in those same units and the output stochastic water paths as well.

      real,    intent(in) :: ciwp    (nlay,dncol)  ! In-cloud ice water path
      real,    intent(in) :: clwp    (nlay,dncol)  ! In-cloud liquid water path

      ! Tiny threshold value for cloud water path, at or below which each of
      ! ice and liquid water parths are separately reset to zero. If both are
      ! thus reset, then the subcolumn gridcell is reset to .not.cloudy. The
      ! idea here is mainly efficiency, since external processing of cloudy
      ! cells usually takes longer.

      real,    intent(in) :: cwp_tiny

      ! output subcolumns
      ! (units of water paths are the same as for inputs ciwp and clwp)
      logical, intent(out) :: cldy_stoch (nlay,nsubcol,dncol)  ! Cloudy or not?
      real,    intent(out) :: ciwp_stoch (nlay,nsubcol,dncol)  ! In-cloud ice water path
      real,    intent(out) :: clwp_stoch (nlay,nsubcol,dncol)  ! In-cloud liq water path

      ! optional seed_order as described above in header
      integer, intent(in), optional :: seed_order (4)

      ! ----- Locals -----

      ! correlation length for cloud presence and condensate
      real, dimension(dncol) :: adl, rdl

      ! inter-layer correlations for cloud presence and condensate
      real, dimension(nlay) :: alpha, rcorr

      ! condensate inhomogeneity
      real :: sigma_qcw(nlay)

      ! seeds for rng_kiss
      real :: pseed(4)  ! isolate seeding pressures
      integer :: so(4)  ! actual seed order used
      logical :: hit(4) ! used in so validation
      integer :: seed1, seed2, seed3, seed4
      integer, parameter :: maximo = huge(seed1) - 1

      ! random number arrays used for overlap
      real, dimension(nlay) :: &
         cdf1, &  ! for cloud presence
         cdf2, &  ! auxilliary
         cdf3     ! for cloud condensate

      ! other locals
      real :: zcw
      logical :: cond_inhomo, ciwp_negligible, clwp_negligible

      ! indices
      integer :: icol, isubcol, ilay, n

      ! save for speed
      cond_inhomo = condensate_inhomogeneous()

      ! set seed order to use
      if (present(seed_order)) then
        do n = 1,4
          so(n) = seed_order(n)
          hit(n) = .false.
        end do
        ! validate a permutation of [1,2,3,4]
        do n = 1,4
          if (so(n) < 1) then
            write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
            error stop 'seed_order element < 1'
          end if
          if (so(n) > 4) then
            write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
            error stop 'seed_order element > 4'
          end if
          if (.not.hit(n)) then
            hit(n) = .true.
          else
            write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
            error stop 'seed_order repeated element'
          end if
        end do
        ! now have 4 elements 1--4 with no repeats so we are good!
      else
        do n = 1,4
          so(n) = n
        end do
      end if

      ! Compute decorrelation length scales ...

      ! for cloud presence
      call correlation_length(dncol, ncol, &
         aam1, aam2, aam30, aam4, doy, alat, adl)

      ! for cloud condensate
      if (cond_inhomo) then
         call correlation_length(dncol, ncol, &
            ram1, ram2, ram30, ram4, doy, alat, rdl)
      endif

      ! outer loop over gridcolumn
      do icol = 1,ncol

         ! exponential inter-layer correlations ...
         do ilay = 2,nlay
            alpha(ilay) = exp( -(zmid(ilay,icol)-zmid(ilay-1,icol)) / adl(icol) )
         enddo
         if (cond_inhomo) then
            do ilay = 2,nlay
               rcorr(ilay) = exp( -(zmid(ilay,icol)-zmid(ilay-1,icol)) / rdl(icol) )
            enddo
         endif

         ! precalculate level of condensate inhomogeneity based on cldfrac.
         ! put outside of subcol loop for speed even though only needed for
         ! cloudy subcolumn cells. sigma_qcw is dimensioned nlay.
         if (cond_inhomo) then
            where (cldfrac(:,icol) > 0.99)
               sigma_qcw = 0.5
            elsewhere (cldfrac(:,icol) > 0.9)
               sigma_qcw = 0.71
            elsewhere
               sigma_qcw = 1.0
            endwhere
         endif

         ! Choose four i4 (32-bit) seeds for the KISS PRNG ...
         !   (a) To ensure reproducibility, choose seeds based on the model state,
         ! in this case based on gridcolumn near surface pressure (this is to be
         ! contrasted with seeds based on, e.g., the system clock, which will give
         ! non-reproducible runs). This state-based reproducibility will regress
         ! with changes to the domain decomposition among processors (so long as
         ! the rest of the model regresses to such changes) or with changes to
         ! the RRTMG gridcolumn partitioning, etc.
         !   (b) If two gridcolumns (or one gridcolumn at two times) have exactly
         ! the same seed, they will produce the same KISS pseudo-random number
         ! sequence. But even if the seeds are different at all, they will still
         ! produce different random sequences. The geostrophic pressure gradient
         ! is ~ rho f u, where f is the coriolis parameter, about 10^-4 s^-1 at
         ! mid-latitudes and u is the wind speed along the isobars. So for u ~ 10
         ! m/s, the pressure gradient is about 1 Pa / km. Of course this varies
         ! with latitude and will not apply in the Tropics (where geostrophy does
         ! not hold). Still it gives us a rough estimate. Now for the high-resol-
         ! ution runs beginning to run, the grid spacing is approaching 1 km, or
         ! about 1 Pa difference! And if two neighboring gridcolumns fall along
         ! an isobar, the pressure difference could be much less. This is why we
         ! would not use the integer portion of the surface pressure in Pa. But
         ! the fractional part of this pressure (scaled to a 32-bit integer) will
         ! yield different seeds.
         !   (c) Having seeds local to each gridcolumn allows potential parallel-
         ! ism. And because the seed is based on pressure differences, not cloud
         ! differences, two gridcolumns adjacent in space or time that have very
         ! similar cloud properties will still generate different subcolumn cloud
         ! ensembles because of their different seeds. This will help beat down
         ! the sampling errors due to finite nsubcol (i.e., subcolumn / g-point
         ! number) when averaging over time and/or space.
         !   (d) Concerning the pressures for the lowest four layers used to
         ! produce the four required KISS seeds: dp ~ rho g dz ~ 10 dz, so a 1 Pa
         ! difference occurs in only 10 cm, which is way smaller than our layer
         ! spacing. In fact, the lowest layers are spaced by about 15 hPa or 1500
         ! Pa for 72L and 2-3 hPa or 200-300 Pa for 181L, both well in excess of
         ! 1 Pa. So we dont need to worry about vertical correlations among the
         ! pseed.

         ! isolate lower pressures for seeding
         pseed = play(1:4,icol) * 100.  ! [Pa]

         ! PMN 2021-11-12: Try [daPa=decaPascals] instead.
         ! pseed = play(1:4,icol) * 10.  ! [daPa]
         ! I was potentially concerned that the fractional part of pseed in [Pa]
         ! would have insufficient significant digits since <play> is 32-bit real.
         ! So, for e.g., 99213.xx [Pa], since 32-bit reals have ~ 7 significant
         ! digits. This would be 9921.3xx [daPa]. But a test of this change made
         ! a minimal change ~0.01 W/m2 in global mean TOA and SFC fluxes and was
         ! not visible in absolute zonal plots. The mapped differences were of
         ! the same nature as other minor changes, perhaps slightly worse. So
         ! I decided not keep this change.

         ! Scaling to integer seeds ...
         ! The 32-bit integer range is [-2147483648,2147483647] and we wish
         ! to avoid zero seeds, so we scale the seeds to [1,2147483647] with
         ! the following (made more general with maximo = huge(integer)-1):
         seed1 = (pseed(so(1)) - int(pseed(so(1)))) * maximo + 1
         seed2 = (pseed(so(2)) - int(pseed(so(2)))) * maximo + 1
         seed3 = (pseed(so(3)) - int(pseed(so(3)))) * maximo + 1
         seed4 = (pseed(so(4)) - int(pseed(so(4)))) * maximo + 1

         ! Generate each subcolumn ...
         do isubcol = 1,nsubcol

            ! exponential overlap in cloud presence
            do ilay = 1,nlay
               call rng_kiss(seed1,seed2,seed3,seed4,cdf1(ilay))
               call rng_kiss(seed1,seed2,seed3,seed4,cdf2(ilay))
            enddo
            do ilay = 2,nlay
               if (cdf2(ilay) < alpha(ilay)) then
                  cdf1(ilay) = cdf1(ilay-1)
               endif
            enddo

            if (cond_inhomo) then

               ! exponential overlap in condensate
               do ilay = 1,nlay
                  call rng_kiss(seed1,seed2,seed3,seed4,cdf2(ilay))
                  call rng_kiss(seed1,seed2,seed3,seed4,cdf3(ilay))
               enddo
               do ilay = 2,nlay
                  if (cdf2(ilay) < rcorr(ilay)) then
                     cdf3(ilay) = cdf3(ilay-1)
                  endif
               enddo

            endif

            ! generate layers of subcolumn ...

            do ilay = 1,nlay

               if (cdf1(ilay) >= 1. - cldfrac(ilay,icol)) then

                  ! a cloudy subcolumn/layer
                  ! note: cdf1 from random num in (0,1) so:
                  !       cldfrac == 0. can never get here;
                  !       cldfrac == 1. always comes here.

                  cldy_stoch(ilay,isubcol,icol) = .true.

                  if (cond_inhomo) then

                     ! horizontally variable clouds
                     zcw = zcw_lookup(cdf3(ilay),sigma_qcw(ilay))
                     ciwp_stoch(ilay,isubcol,icol) = ciwp(ilay,icol) * zcw
                     clwp_stoch(ilay,isubcol,icol) = clwp(ilay,icol) * zcw

                  else

                     ! horizontally homogeneous clouds
                     ciwp_stoch(ilay,isubcol,icol) = ciwp(ilay,icol)
                     clwp_stoch(ilay,isubcol,icol) = clwp(ilay,icol)

                  endif

                  ! reset negligible water paths to zero (see comment in intro)
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

                  cldy_stoch(ilay,isubcol,icol) = .false.
                  ciwp_stoch(ilay,isubcol,icol) = 0.
                  clwp_stoch(ilay,isubcol,icol) = 0.

               endif
           
            enddo  ! layer

         enddo  ! subcolumn
      enddo  ! gridcolumn

   end subroutine generate_stochastic_clouds


   !-------------------------------------------
   subroutine correlation_length(dncol, ncol, &
      am1, am2, am30, am4, doy, alat, clength)
   !-------------------------------------------
      integer, intent(in)  :: dncol                ! Dimensioned number of gridcols
      integer, intent(in)  :: ncol                 ! Actual number of gridcols
      real,    intent(in)  :: am1, am2, am30, am4  ! input parameters
      integer, intent(in)  :: doy                  ! Day of year
      real,    intent(in)  :: alat    (dncol)      ! Latitude of gridcolumn [radians]
      real,    intent(out) :: clength (dncol)      ! Correlation length [m]

      real, parameter :: r2d = 180.d0 / 3.14159265358979323846d0

      integer :: icol
      real :: am3

      if (doy > 181) then
         am3 = -4.*am30/365.*(doy-272)
      else
         am3 =  4.*am30/365.*(doy- 91)
      endif

      do icol = 1,ncol
         clength(icol) = (am1+am2*exp(-(alat(icol)*r2d-am3)**2/am4**2))*1.e3
      enddo

   end subroutine correlation_length


   !-------------------------------------------------------
   subroutine rng_kiss(seed1, seed2, seed3, seed4, ran_num)
   !-------------------------------------------------------
   ! See note below: get ran_num in (0,1).
   !
   ! public domain code.
   ! made available from http://www.fortran.com/.
   ! downloaded by pjr on 03/16/04 for NCAR CAM.
   ! converted to vector form, functions inlined by pjr,mvr on 05/10/2004.
   ! devectorized by pmn 08/24/2021.
   !
   ! The KISS (Keep It Simple Stupid) random number generator. Combines:
   ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
   ! (2) A 3-shift shift-register generator, period 2^32-1,
   ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
   ! Overall period>2^123.
   !-------------------------------------------------------

      integer, intent(inout) :: seed1, seed2, seed3, seed4
      real, intent(out) :: ran_num
      integer :: m, k, n, kiss

      ! inline function 
      m(k,n) = ieor(k,ishft(k,n))

      seed1 = 69069 * seed1 + 1327217885
      seed2 = m(m(m(seed2,13),-17),5)
      seed3 = 18000 * iand(seed3,65535) + ishft(seed3,-16)
      seed4 = 30903 * iand(seed4,65535) + ishft(seed4,-16)
      kiss = seed1 + seed2 + ishft(seed3,16) + seed4
      ran_num = kiss * 2.328306e-10 + 0.5 
    
   !-------------------------------------------------------
   ! pmn notes:
   ! kiss is a random 32-bit signed integer in [-(2**31),2**31-1],
   ! and 2.328306e-10 is an approx to 2**-32 =~ 2.32830644e-10.
   ! Based on the tests below, the result ran_num will be
   !    [8.9406967E-08,0.9999999] ~ (0,1).
   !
   ! rng_test.f90:
   ! program rng_test
   !    integer, parameter :: mini = -2**31
   !    integer, parameter :: maxi = 2**31 - 1
   !    real :: rng
   !    write(*,*) mini
   !    write(*,*) maxi
   !    rng = mini * 2.328306e-10 + 0.5
   !    write(*,*) rng
   !    rng = maxi * 2.328306e-10 + 0.5
   !    write(*,*) rng
   ! end program rng_test
   !
   ! Results:
   !  $ module load comp/intel/19.1.3.304
   !  $ ifort rng_test.f90
   !  $ a.out
   !  -2147483648
   !   2147483647
   !   8.9406967E-08
   !   0.9999999
   !-------------------------------------------------------

   end subroutine rng_kiss


   ! ----------------------------------------------------------------
   subroutine clearCounts_threeBand( &
                 dncol, ncol, nsubcol, nlay, &
                 cloudLM, cloudMH, cldy_stoch, &
                 clearCnts)
   ! ----------------------------------------------------------------
   ! Count clear subcolumns per gridcolumn for whole column
   !   and for three pressure bands.
   ! layers [1,         cloudLM] are in low  pressure band
   ! layers [cloudLM+1, cloudMH] are in mid  pressure band
   ! layers [cloudMH+1, nlay   ] are in high pressure band
   ! ----------------------------------------------------------------

      integer, intent(in)  :: dncol              ! Dimensioned number of gridcols
      integer, intent(in)  :: ncol               ! Actual number of gridcols
      integer, intent(in)  :: nsubcol            ! number of subcols per gridcol
      integer, intent(in)  :: nlay               ! number of layers
      integer, intent(in)  :: cloudLM, cloudMH   ! layer indices as above

      logical, intent(in)  :: cldy_stoch(nlay,nsubcol,dncol)  ! cloudy or not?

      integer, intent(out) :: clearCnts(4,dncol) ! counts of clear subcolumns
                                                 ! (1,) whole column
                                                 ! (2,) high pressure band only
                                                 ! (3,) mid  pressure band only
                                                 ! (4,) low  pressure band only
                                          
      ! locals
      integer :: icol, isubcol, ilay
      logical :: cloud_found

      ! zero counters
      clearCnts = 0
 
      do icol = 1,ncol
         do isubcol = 1,nsubcol

            ! whole subcolumn
            cloud_found = .false.
            do ilay = 1, nlay
               if (cldy_stoch(ilay,isubcol,icol)) then 
                  cloud_found = .true.
                  exit
               endif
            enddo
            if (.not. cloud_found) &
               clearCnts(1,icol) = clearCnts(1,icol) + 1

            ! high pressure band
            cloud_found = .false.
            do ilay = cloudMH+1, nlay
               if (cldy_stoch(ilay,isubcol,icol)) then 
                  cloud_found = .true.
                  exit
               endif
            enddo
            if (.not. cloud_found) &
               clearCnts(2,icol) = clearCnts(2,icol) + 1

            ! mid pressure band
            cloud_found = .false.
            do ilay = cloudLM+1, cloudMH
               if (cldy_stoch(ilay,isubcol,icol)) then 
                  cloud_found = .true.
                  exit
               endif
            enddo
            if (.not. cloud_found) &
               clearCnts(3,icol) = clearCnts(3,icol) + 1

            ! low pressure band
            cloud_found = .false.
            do ilay = 1, cloudLM
               if (cldy_stoch(ilay,isubcol,icol)) then 
                  cloud_found = .true.
                  exit
               endif
            enddo
            if (.not. cloud_found) &
               clearCnts(4,icol) = clearCnts(4,icol) + 1

         enddo  ! subcolumns
      enddo  ! gridcolumns
   
   end subroutine clearCounts_threeBand

end module cloud_subcol_gen
