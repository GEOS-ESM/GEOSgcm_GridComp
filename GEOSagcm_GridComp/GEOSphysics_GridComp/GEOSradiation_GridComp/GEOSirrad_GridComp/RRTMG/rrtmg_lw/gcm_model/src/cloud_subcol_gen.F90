module cloud_subcol_gen

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! Purpose: Create stochastic arrays for cloud physical properties.
! Input gridcolumn cloud profiles: cloud fraction and in-cloud ice and liquid
! water paths. Output will be stochastic subcolumn arrays of these variables.

! --------- Module variables ----------

      use cloud_condensate_inhomogeneity, only: condensate_inhomogeneous, zcw_lookup

      implicit none
      private

      public :: generate_stochastic_clouds 
      public :: clearCounts_threeBand

contains

!----------------------------------------------------------------------------------------------
   subroutine generate_stochastic_clouds( &
                 ncol, nsubcol, nlay, &
                 zmid, alat, doy, &
                 play, cldfrac, ciwp, clwp, &
                 cldf_stoch, ciwp_stoch, clwp_stoch)
!----------------------------------------------------------------------------------------------
!
! Original code: Based on Raisanen et al., QJRMS, 2004.
! Original Contact: Cecile Hannay (hannay@ucar.edu)
! 
! Modifications:
! ...
! Peter Norris, GMAO, Apr 2021:
!   Tidy up code and comments
!   Keep only exponential (generalized) overlap.
!
! Given a profile of cloud fraction, cloud water and cloud ice, produce a set of subcolumns.
! Each subcolumn has cloud fraction in {0,1} at each level.
! If (.not. condensate_inhomogeneous()), each layer within each subcolumn has uniform
!   cloud ice and cloud liquid concentration.
! If (condensate_inhomogeneous()), each layer has horizontal condensate variability.
! The ensemble of subcolumns statistically reproduces the cloud fraction within each layer
!   and its prescribed vertical overlap, and, if (condensate_inhomogeneous()), the PDF of
!   cloud liquid and ice within each layer and its prescribed vertical correlation.
! Whether there is condensate inhomogeneity, and the condensate PDF if so, are set in
!   initialize_inhomogeneity().
! 
! Overlap assumption:
! Exponential (generalized) overlap (Raisanen et al., 2004, Pincus et al., 2005) using
! correlations alpha and rcorr based on decorrelation length scales from Oreopoulos et
! al., 2012.
! 
!---------------------------------------------------------------------------------------------

      integer, intent(in) :: ncol                 ! number of gridcolumns
      integer, intent(in) :: nsubcol              ! number of subcols to generate / gridcol
      integer, intent(in) :: nlay                 ! number of model layers
      real,    intent(in) :: zmid    (ncol,nlay)  ! Hgt of midpoints above sfc [m]
      real,    intent(in) :: alat    (ncol)       ! Latitude of gridcolumn
      integer, intent(in) :: doy                  ! Day of year
      real,    intent(in) :: play    (ncol,nlay)  ! layer pressures [mb]
      real,    intent(in) :: cldfrac (ncol,nlay)  ! layer cloud fraction
      real,    intent(in) :: ciwp    (ncol,nlay)  ! in-cld ice water path
      real,    intent(in) :: clwp    (ncol,nlay)  ! in-cld liquid water path

      ! output subcolumns, one subcolum per gpoint
      real,    intent(out) :: cldf_stoch(ncol,nsubcol,nlay)  ! cloud fraction 
      real,    intent(out) :: ciwp_stoch(ncol,nsubcol,nlay)  ! in-cloud ice water path
      real,    intent(out) :: clwp_stoch(ncol,nsubcol,nlay)  ! in-cloud liq water path
   
      ! ----- Locals -----

      ! layer pressures [Pa] used for seeding
      real :: pmid(ncol,nlay)

      ! decorrelation length scales for cldfrac and condensate
      real :: adl(ncol), rdl(ncol)

      ! inter-layer correlations for cldfrac and condensate
      real :: alpha(nlay), rcorr(nlay)

      ! related to decorrelation lengths
      real :: am1, am2, am3, am4, amr
  
      ! related to random number and seed 
      integer :: seed1, seed2, seed3, seed4   ! seed (kissvec)
      real :: rand_num                        ! random number (kissvec)

      ! random number arrays used for overlap
      real :: cdf1(nlay)  ! for cloud presence
      real :: cdf2(nlay)  ! auxilliary
      real :: cdf3(nlay)  ! for cloud condensate

      ! other locals
      real :: cfs, zcw, sigma_qcw
      logical :: cond_inhomo
     
      ! indices
      integer :: icol, isubcol, ilev
      
      ! midpoint pressure in [Pa] for seeding
      pmid = play * 1.e2

      ! save for speed
      cond_inhomo = condensate_inhomogeneous()

      ! -----------------------------------
      ! compute decorrelation length scales
      ! -----------------------------------

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
      adl = (am1+am2*exp(-(alat*180./3.141592-am3)**2/am4**2))*1.e3  ! [m]

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
      rdl = (am1+am2*exp(-(alat*180./3.141592-am3)**2/am4**2))*1.e3  ! [m]
        
      ! --------------------------
      ! outer loop over gridcolumn
      ! --------------------------
      do icol = 1, ncol

         ! ------------------------------------
         ! exponential inter-layer correlations
         ! ------------------------------------
         alpha(1) = 0.0
         rcorr(1) = 0.0
         do ilev = 2, nlay
            alpha(ilev) = exp( -(zmid(icol,ilev)-zmid(icol,ilev-1)) / adl(icol) )
            rcorr(ilev) = exp( -(zmid(icol,ilev)-zmid(icol,ilev-1)) / rdl(icol) )
         end do

         ! -----------------------
         ! generate each subcolumn
         ! -----------------------
         do isubcol = 1, nsubcol

            ! -----------
            ! create seed
            ! -----------
            ! For kissvec, create a seed that depends on the state of the columns.
            ! Maybe not the best way, but it works. Must use pmid from bottom four layers. 
            ! pmn ... why do we have to reset seed at beginning of every subcol?
            ! pmn ... no real harm and allowed parallel generation in old gpu code
      
            seed1 = (pmid(icol,1) - int(pmid(icol,1))) * 1000000000 + isubcol * 11 
            seed3 = (pmid(icol,3) - int(pmid(icol,3))) * 1000000000 + isubcol * 13 
            seed2 = seed1 + isubcol
            seed4 = seed3 - isubcol 
  
            ! ------------------------
            ! apply overlap assumption
            ! ------------------------

            do ilev = 1,nlay
               call kissvec(seed1,seed2,seed3,seed4,rand_num)
               cdf1(ilev) = rand_num
               call kissvec(seed1,seed2,seed3,seed4,rand_num)
               cdf2(ilev) = rand_num
            end do

            ! exponential overlap in cloud fraction
            do ilev = 2,nlay
               if (cdf2(ilev) < alpha(ilev)) then
                  cdf1(ilev) = cdf1(ilev-1)
               end if
            end do

            ! exponential overlap in condensate
            if (cond_inhomo) then

               do ilev = 1,nlay
                  call kissvec(seed1,seed2,seed3,seed4,rand_num)
                  cdf2(ilev) = rand_num
                  call kissvec(seed1,seed2,seed3,seed4,rand_num)
                  cdf3(ilev) = rand_num 
               end do

               do ilev = 2,nlay
                  if (cdf2(ilev) < rcorr(ilev)) then
                     cdf3(ilev) = cdf3(ilev-1) 
                  end if
               end do

            end if

            ! -------------------
            ! generate subcolumns
            ! -------------------

            do ilev = 1,nlay
               cfs = cldfrac(icol,ilev)

               ! if a cloudy subcolumn
               if (cdf1(ilev) >= 1. - cfs) then

                  cldf_stoch(icol,isubcol,ilev) = 1. 

                  ! condensate is horizontally homogeneous by default
                  clwp_stoch(icol,isubcol,ilev) = clwp(icol,ilev)
                  ciwp_stoch(icol,isubcol,ilev) = ciwp(icol,ilev)
          
                  ! horizontal condensate variability if requested
                  if (cond_inhomo) then  

                     ! Cloud fraction sets level of inhomogeneity
                     if (cfs .gt. 0.99 ) then
                        sigma_qcw = 0.5
                     elseif (cfs .gt. 0.9 ) then
                        sigma_qcw = 0.71
                     else  
                        sigma_qcw = 1.0
                     endif

                     ! horizontally variable clouds
                     zcw = zcw_lookup(cdf3(ilev),sigma_qcw)
                     clwp_stoch(icol,isubcol,ilev) = clwp(icol,ilev) * zcw
                     ciwp_stoch(icol,isubcol,ilev) = ciwp(icol,ilev) * zcw
                  
                  end if

               else

                  ! a clear subcolumn
                  cldf_stoch(icol,isubcol,ilev) = 0. 
                  clwp_stoch(icol,isubcol,ilev) = 0. 
                  ciwp_stoch(icol,isubcol,ilev) = 0. 

               endif
           
            end do  ! level

         end do  ! isubcol
      end do  ! icol

   end subroutine generate_stochastic_clouds


!------------------------------------------------------------
   subroutine kissvec(seed1, seed2, seed3, seed4, ran_arr)
!------------------------------------------------------------
!
! public domain code
! made available from http://www.fortran.com/
! downloaded by pjr on 03/16/04 for NCAR CAM
! converted to vector form, functions inlined by pjr,mvr on 05/10/2004
!
! The KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
! Overall period>2^123; 
!------------------------------------------------------------

      real, intent(inout) :: ran_arr
      integer, intent(inout) :: seed1, seed2, seed3, seed4
      integer :: m, k, n, kiss

      ! inline function 
      m(k, n) = ieor (k, ishft (k, n) )

      seed1 = 69069 * seed1 + 1327217885
      seed2 = m (m (m (seed2, 13), - 17), 5)
      seed3 = 18000 * iand (seed3, 65535) + ishft (seed3, - 16)
      seed4 = 30903 * iand (seed4, 65535) + ishft (seed4, - 16)
      kiss = seed1 + seed2 + ishft (seed3, 16) + seed4
      ran_arr = kiss*2.328306e-10 + 0.5 
    
   end subroutine kissvec


! ------------------------------------------------------------------------------------
   subroutine clearCounts_threeBand( &
                 ncol, nsubcol, nlay, cloudLM, cloudMH, cldf_stoch, &
                 clearCounts)
! ------------------------------------------------------------------------------------
!
! count clear subcolumns per gridcolumn for whole column and three pressure bands
! layers [1,         cloudLM] are in low  pressure band
! layers [cloudLM+1, cloudMH] are in mid  pressure band
! layers [cloudMH+1, nlay   ] are in high pressure band
! ------------------------------------------------------------------------------------

      integer, intent(in)  :: ncol               ! number of gridcolumns
      integer, intent(in)  :: nsubcol            ! number of subcols per gridcol
      integer, intent(in)  :: nlay               ! number of layers
      integer, intent(in)  :: cloudLM, cloudMH   ! layer indices as above

      real,    intent(in)  :: cldf_stoch(ncol,nsubcol,nlay)  ! cloud fraction [mcica]

      integer, intent(out) :: clearCounts(ncol,4)  ! counts of clear subcolumns
                                                   ! (,1) whole column
                                                   ! (,2) high pressure band only
                                                   ! (,3) mid  pressure band only
                                                   ! (,4) low  pressure band only
                                          
      ! locals
      integer :: icol, isubcol, ilev, tk

      ! zero counters
      clearCounts = 0
 
      do icol = 1, ncol
         do isubcol = 1, nsubcol

            ! whole subcolumn
            ! remember cldf_stoch in {0,1}
            tk = 0
            do ilev = 1, nlay
               if (cldf_stoch(icol,isubcol,ilev) > 0.5) then 
                  tk = 1
                  exit
               end if
            end do
            if (tk .eq. 0) & 
               clearCounts(icol,1) = clearCounts(icol,1) + 1

            ! high pressure band
            tk = 0
            do ilev = cloudMH+1, nlay
               if (cldf_stoch(icol,isubcol,ilev) > 0.5) then 
                  tk = 1
                  exit
               end if
            end do
            if (tk .eq. 0) & 
               clearCounts(icol,2) = clearCounts(icol,2) + 1

            ! mid pressure band
            tk = 0
            do ilev = cloudLM+1, cloudMH
               if (cldf_stoch(icol,isubcol,ilev) > 0.5) then 
                  tk = 1
                  exit
               end if
            end do
            if (tk .eq. 0) & 
               clearCounts(icol,3) = clearCounts(icol,3) + 1

            ! low pressure band
            tk = 0
            do ilev = 1, cloudLM
               if (cldf_stoch(icol,isubcol,ilev) > 0.5) then 
                  tk = 1
                  exit
               end if
            end do
            if (tk .eq. 0) & 
               clearCounts(icol,4) = clearCounts(icol,4) + 1

         end do  ! isubcol
      end do  ! icol
   
   end subroutine clearCounts_threeBand

end module cloud_subcol_gen
