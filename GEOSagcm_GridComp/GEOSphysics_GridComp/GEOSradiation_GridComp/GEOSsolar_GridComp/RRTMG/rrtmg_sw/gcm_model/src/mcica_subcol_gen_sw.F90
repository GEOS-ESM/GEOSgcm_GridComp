module mcica_subcol_gen_sw

   ! Purpose: Create stochastic arrays for cloud physical properties.
   ! Input gridcolumn cloud profiles: cloud fraction and in-cloud ice and liquid
   ! water paths. Output will be stochastic subcolumn arrays of these variables.

   use cloud_condensate_inhomogeneity, only: &
      condensate_inhomogeneous, zcw_lookup
    
   implicit none
   private

   public :: mcica_sw      
      
contains

   !---------------------------------------------------------------------------------------
   subroutine mcica_sw( &
      pncol, ncol, nsubcol, nlay, &
      zmid, alat, doy, &
      play, cld, ciwp, clwp, &
      cld_stoch, ciwp_stoch, clwp_stoch)
   !---------------------------------------------------------------------------------------
   !
   ! Original code: Based on Raisanen et al., QJRMS, 2004.
   ! Original Contact: Cecile Hannay (hannay@ucar.edu)
   !
   ! Modifications:
   ! ...
   ! Peter Norris, GMAO, Apr-May 2021:
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

      integer, intent(in) :: pncol              ! Dimensioned number of gridcols
      integer, intent(in) :: ncol               ! Actual number of gridcols
      integer, intent(in) :: nsubcol            ! #Subcols to generate / gridcol
      integer, intent(in) :: nlay               ! Number of model layers
      real,    intent(in) :: zmid (nlay,pncol)  ! Hgt of midpoints [m]
      real,    intent(in) :: alat      (pncol)  ! Latitude of gridcolumn
      integer, intent(in) :: doy                ! Day of year
      real,    intent(in) :: play (nlay,pncol)  ! Layer pressures [Pa]
!pmn: these seem to be passed in as hPa!!!!!
      real,    intent(in) :: cld  (nlay,pncol)  ! Layer cloud fraction 
      real,    intent(in) :: ciwp (nlay,pncol)  ! In-cloud ice water path (g/m2)?
      real,    intent(in) :: clwp (nlay,pncol)  ! In-cloud liquid water path (g/m2)?

      ! output subcolumns
      ! (units of water paths are the same as for inputs ciwp and clwp)
      real,    intent(out) ::  cld_stoch (nlay,nsubcol,pncol)  ! Cloud fraction 
      real,    intent(out) :: ciwp_stoch (nlay,nsubcol,pncol)  ! In-cloud ice water path
      real,    intent(out) :: clwp_stoch (nlay,nsubcol,pncol)  ! In-cloud liq water path
      
      ! ----- Locals -----

!pmn?  units and use cf LW
      ! Constants (min value for cloud fraction and cloud water and ice)
      real, parameter :: cldmin = 1.0e-20    ! min cloud fraction

      ! decorrelation length scales for cldfrac and condensate
      real, dimension(pncol) :: adl, rdl

      ! inter-layer correlations for cldfrac and condensate
      real, dimension(pncol,nlay) :: alpha, rcorr

      ! related to decorrelation lengths
      real :: am1, am2, am3, am4, amr

      ! related to random number and seed 
      integer :: seed1, seed2, seed3, seed4  ! seed (rng_kiss)
      real :: rand_num                       ! random number (rng_kiss)

      ! random number arrays used for overlap
      real, dimension(pncol,nlay,nsubcol) :: &
         cdf1, &  ! for cloud presence
         cdf2, &  ! auxilliary
         cdf3     ! for cloud condensate

      ! other locals
      real :: sigma_qcw, zcw !? cfs
      logical :: cond_inhomo !?, ciwp_negligible, clwp_negligible

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

!pmn y
      ! cloud presence decorrelation length scale
      am1 = 1.4315
      am2 = 2.1219
      am4 = -25.584
      amr = 7.
      if (doy .gt. 181) then
         am3 = -4.*amr/365*(doy-272)
      else
         am3 = 4.*amr/365.*(doy-91)
      endif
      do icol = 1,ncol
         adl(icol) = (am1+am2*exp(-(alat(icol)*180./3.141592-am3)**2/am4**2))*1.e3  ! [m]
      end do

!pmn y
      ! condensate decorrelation length scale
      am1 = 0.72192
      am2 = 0.78996
      am4 = 40.404
      amr = 8.5
      if (doy .gt. 181) then
         am3 = -4.*amr/365*(doy-272)
      else
         am3 = 4.*amr/365.*(doy-91)
      endif
      do icol = 1,ncol
         rdl(icol) = (am1+am2*exp(-(alat(icol)*180./3.141592-am3)**2/am4**2))*1.e3  ! [m]
      end do
   
!pmn y
      do icol = 1,ncol
         alpha(icol,1) = 0.
         do ilay = 2,nlay
            alpha(icol,ilay) = exp(-(zmid(ilay,icol) - zmid(ilay-1,icol)) / adl(icol))
         end do
      end do
     
! pmn icol ily nsub ie (nsub,ily,icl)
! but non-zero-diff <<<<<<<<<<<<<<<<<<<<<  !!!!
! preserve zero-diff as longg as possible so work on other parts of code first

      do ilay = 1,nlay
         do icol = 1,ncol
            seed1 = (play(1,icol) - int(play(1,icol))) * 100000000 - ilay
            seed2 = (play(2,icol) - int(play(2,icol))) * 100000000 + ilay
            seed3 = (play(3,icol) - int(play(3,icol))) * 100000000 + ilay * 6.2
            seed4 = (play(4,icol) - int(play(4,icol))) * 100000000           
!pmn: 8 zeros here cf 9 for LW!!!
            do isubcol = 1,nsubcol
               call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
               cdf1(icol,ilay,isubcol) = rand_num 
               call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
               cdf2(icol,ilay,isubcol) = rand_num 
            end do
         end do
      end do
    
!pmn y (isub,ilay,icol)
      do icol = 1,ncol
         do isubcol = 1,nsubcol
            do ilay = 2,nlay
               if (cdf2(icol,ilay,isubcol) < alpha(icol,ilay)) then
                  cdf1(icol,ilay,isubcol) = cdf1(icol,ilay-1,isubcol) 
               end if
            end do
         end do
      end do
    
      if (cond_inhomo) then
        
         do icol = 1,ncol
            rcorr(icol,1) = 0.
            do ilay = 2,nlay
               rcorr(icol,ilay) = exp(-(zmid(ilay,icol) - zmid(ilay-1,icol)) / rdl(icol))
            end do
         end do
        
         do ilay = 1,nlay
            do icol = 1,ncol
               seed1 = (play(1,icol) - int(play(1,icol))) * 100000000 - ilay
               seed2 = (play(2,icol) - int(play(2,icol))) * 100000000 + ilay
               seed3 = (play(3,icol) - int(play(3,icol))) * 100000000 + ilay * 6.2
               seed4 = (play(4,icol) - int(play(4,icol))) * 100000000           
               do isubcol = 1,nsubcol
                  call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
                  cdf2(icol,ilay,isubcol) = rand_num 
                  call rng_kiss(seed1,seed2,seed3,seed4,rand_num)
                  cdf3(icol,ilay,isubcol) = rand_num 
               end do
            end do
         end do

         do icol = 1,ncol
            do isubcol = 1,nsubcol
               do ilay = 2,nlay
                  if (cdf2(icol,ilay,isubcol) < rcorr(icol,ilay)) then
                     cdf3(icol,ilay,isubcol) = cdf3(icol,ilay-1,isubcol)
                  end if
               end do
            end do
         end do
      end if

      ! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
      ! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0;
      ! where there is a cloud, define the subcolumn cloud properties,
      ! otherwise set these to zero

!pmn refactor of cond_inhomo
! otherwise I think y
      do ilay = 1,nlay
         do icol = 1,ncol
            do isubcol = 1,nsubcol

               if (cond_inhomo .and. cdf1(icol,ilay,isubcol) >= (1. - cld(ilay,icol))) then

                  cld_stoch(ilay,isubcol,icol) = 1. 
                  
                  ! Cloud fraction sets level of inhomogeneity
                  if (cld(ilay,icol) .gt. 0.99) then
                     sigma_qcw = 0.5
                  elseif (cld(ilay,icol) .gt. 0.9) then
                     sigma_qcw = 0.71
                  else  
                     sigma_qcw = 1.0
                  endif
                  
                  ! horizontally variable clouds
                  zcw = zcw_lookup(cdf3(icol,ilay,isubcol),sigma_qcw)
                  clwp_stoch(ilay,isubcol,icol) = clwp(ilay,icol) * zcw
                  ciwp_stoch(ilay,isubcol,icol) = ciwp(ilay,icol) * zcw
                
               elseif (cdf1(icol,ilay,isubcol) >= (1. - cld(ilay,icol))) then

                   cld_stoch(ilay,isubcol,icol) = 1. 
                  clwp_stoch(ilay,isubcol,icol) = clwp(ilay,icol)
                  ciwp_stoch(ilay,isubcol,icol) = ciwp(ilay,icol)

               else

                  ! a clear subcolumn
                   cld_stoch(ilay,isubcol,icol) = 0. 
                  clwp_stoch(ilay,isubcol,icol) = 0. 
                  ciwp_stoch(ilay,isubcol,icol) = 0. 

               endif
            enddo
         enddo
      enddo

! ilay, icol: alpha, rcorr
! isubcol, ilay, icol: stoch, and perhaps cdf* ?

   end subroutine mcica_sw

   !-------------------------------------------------------
   subroutine rng_kiss(seed1, seed2, seed3, seed4, ran_num)
   !-------------------------------------------------------

      real, intent(inout) :: ran_num
      integer, intent(inout) :: seed1, seed2, seed3, seed4
!     integer :: m, k, n, kiss
      integer :: kiss

!     ! inline function
!     m(k, n) = ieor (k, ishft (k, n) )

     !seed1 = 69069 * seed1 + 1327217885 ! in LW
      seed1 = 69069 * seed1 + 132721785  ! in SW
     !seed2 = m (m (m (seed2, 13), - 17), 5)                    ! in LW
      seed2 = 11002 * iand (seed2, 65535) + ishft (seed2, -16)  ! in SW
      seed3 = 18000 * iand (seed3, 65535) + ishft (seed3, -16)
      seed4 = 30903 * iand (seed4, 65535) + ishft (seed4, -16)
      kiss = seed1 + seed2 + ishft (seed3, 16) + seed4
      ran_num = kiss * 2.328306e-10 + 0.5

   end subroutine rng_kiss

end module mcica_subcol_gen_sw
