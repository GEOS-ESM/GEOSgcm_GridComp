!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
#include "_gpudef.inc"
! (dmb 2012) This is the GPU version of the cldprmc routine.  I have parallelized across 
! all 3 dimensions (columns, g-points, and layers) to make this routine run very fast on the GPU.  
! The greatest speedup was obtained by switching the indices for the cloud variables so that 
! the columns were the least significant (leftmost) dimension
      module gpu_rrtmg_lw_cldprmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! --------- Modules ----------

       !use parkind, only : im => kind , rb => kind 
      use parrrtm, only : ngptlw, nbndlw
      use rrlw_cld, only: abscld1, absliq0, absliq1, &
                          absice0, absice1, absice2, absice3, absice4
     ! use rrlw_wvn, only: ngb
      use rrlw_vsn, only: hvrclc, hnamclc

#ifdef _CUDA
	  use cudafor
#endif
      implicit none

   
	  ! (dmb 2012) I moved most GPU variables so that they are module level variables.
      ! PGI Fortran seems to sometimes have trouble passing arrays into kernels correctly.
      ! Using module level variables bypasses this issue and allows for cleaner code.
      integer  _gpudev, allocatable :: inflagd(:), iceflagd(:), liqflagd(:)

      real  _gpudev, allocatable :: ciwpmcd(:,:,:)  ! in-cloud ice water path [mcica]
      real  _gpudev, allocatable :: clwpmcd(:,:,:)  ! in-cloud liquid water path [mcica]
      real  _gpudev, allocatable :: relqmcd(:,:)         ! liquid particle effective radius (microns)
      real  _gpudev, allocatable :: reicmcd(:,:)         ! ice particle effective size (microns)
 
      real   _gpucon, dimension(2) :: absice0d
      real   _gpucon, dimension(2,5) :: absice1d
      real   _gpucon, dimension(43,16) :: absice2d
      real   _gpucon, dimension(46,16) :: absice3d
      real   _gpucon, dimension(200,16) :: absice4d
      real   _gpucon, dimension(58,16) :: absliq1d


      contains


	  
	  ! ------------------------------------------------------------------------------
      _gpuker subroutine cldprmcg(ncol, nlayers, cldfmc, &
                          taucmc, ngb, icb, ncbands, icldlyr)
! ------------------------------------------------------------------------------

! Purpose:  Compute the cloud optical depth(s) for each cloudy layer.

! ------- Input -------

	  integer, value, intent(in) :: ncol				! total number of columns
      integer, value, intent(in) :: nlayers         ! total number of layers
      
       real , intent(in) :: cldfmc(ncol, ngptlw, nlayers+1)        ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
                                                       ! specific definition of reicmc depends on setting of iceflag:
                                                      ! iceflag = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec must be >= 10.0 microns
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]

	  integer   , intent(out) :: icldlyr( ncol, nlayers+1)
	  integer   , dimension(140), intent(in)  :: ngb
	  integer, intent(in) :: icb(16)
      real , intent(inout) :: taucmc(:,:,:)     ! cloud optical depth [mcica]
	  real , parameter :: absliq0 = 0.0903614 

! ------- Output -------

      integer , intent(out) :: ncbands(:)        ! number of cloud spectral bands
      
                                                      !    Dimensions: (ngptlw,nlayers)

! ------- Local -------

      integer  :: iplon
	  integer  :: lay                         ! Layer index
      integer  :: ib                          ! spectral band index
      integer  :: ig                          ! g-point interval index
      integer  :: index 
     

      real  :: abscoice                       ! ice absorption coefficients
      real  :: abscoliq                       ! liquid absorption coefficients
      real  :: cwp                            ! cloud water path
      real  :: radice                         ! cloud ice effective size (microns)
      real  :: factor                         ! 
      real  :: fint                           ! 
      real  :: radliq                         ! cloud liquid droplet radius (microns)
      real , parameter :: eps = 1.e-6       ! epsilon
      real , parameter :: cldmin = 1.e-20   ! minimum value for cloud quantities

! ------- Definitions -------

!     Explanation of the method for each value of INFLAG.  Values of
!     0 or 1 for INFLAG do not distingish being liquid and ice clouds.
!     INFLAG = 2 does distinguish between liquid and ice clouds, and
!     requires further user input to specify the method to be used to 
!     compute the aborption due to each.
!     INFLAG = 0:  For each cloudy layer, the cloud fraction and (gray)
!                  optical depth are input.  
!     INFLAG = 1:  For each cloudy layer, the cloud fraction and cloud
!                  water path (g/m2) are input.  The (gray) cloud optical 
!                  depth is computed as in CCM2.
!     INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud 
!                  water path (g/m2), and cloud ice fraction are input.
!       ICEFLAG = 0:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in CCM3.
!       ICEFLAG = 1:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in 
!                     Ebert and Curry, JGR, 97, 3831-3836 (1992).  The 
!                     spectral regions in this work have been matched with
!                     the spectral bands in RRTM to as great an extent 
!                     as possible:  
!                     E&C 1      IB = 5      RRTM bands 9-16
!                     E&C 2      IB = 4      RRTM bands 6-8
!                     E&C 3      IB = 3      RRTM bands 3-5
!                     E&C 4      IB = 2      RRTM band 2
!                     E&C 5      IB = 1      RRTM band 1
!       ICEFLAG = 2:  The ice effective radius (microns) is input and the
!                     optical properties due to ice clouds are computed from
!                     the optical properties stored in the RT code,
!                     STREAMER v3.0 (Reference: Key. J., Streamer 
!                     User's Guide, Cooperative Institute for
!                     Meteorological Satellite Studies, 2001, 96 pp.).
!                     Valid range of values for re are between 5.0 and
!                     131.0 micron.
!       ICEFLAG = 3: The ice generalized effective size (dge) is input
!                    and the optical properties, are calculated as in
!                    Q. Fu, J. Climate, (1998). Q. Fu provided high resolution
!                    tables which were appropriately averaged for the
!                    bands in RRTM_LW.  Linear interpolation is used to
!                    get the coefficients from the stored tables.
!                    Valid range of values for dge are between 5.0 and
!                    140.0 micron.
!       LIQFLAG = 0:  The optical depths due to water clouds are computed as
!                     in CCM3.
!       LIQFLAG = 1:  The water droplet effective radius (microns) is input 
!                     and the optical depths due to water clouds are computed 
!                     as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
!                     The values for absorption coefficients appropriate for
!                     the spectral bands in RRTM have been obtained for a 
!                     range of effective radii by an averaging procedure 
!                     based on the work of J. Pinto (private communication).
!                     Linear interpolation is used to get the absorption 
!                     coefficients for the input effective radius.

    ! (dmb 2012) Here insead of looping over the column, layer, and band dimensions,
    ! I compute the index for each dimension from the grid and block layout.  This 
    ! function is called once per each thread, and each thread has a unique combination of 
    ! column, layer, and g-point.  

#ifdef _CUDA
    iplon = (blockidx%x-1) * blockdim%x + threadidx%x
	lay = (blockidx%y-1) * blockdim%y + threadidx%y
    ig = (blockidx%z-1) * blockdim%z + threadidx%z
    ! (dmb 2012) Make sure that the column, layer, and g-points are all within the proper
    ! range.  They can be out of range if we select certain block configurations due to 
    ! optimizations.
    if (iplon<=ncol .and. lay<=nlayers .and. ig<=ngptlw) then
#else
    do iplon = 1, ncol
    do lay= 1, nlayers
    do ig= 1, ngptlw
#endif


    ncbands(iplon) = 1
	      ! (dmb 2012) all of the cloud variables have been modified so that the column dimensions 
          ! is least significant.
          if (cldfmc(iplon,ig,lay) .eq. 1. ) then
            icldlyr(iplon, lay)=1
          endif
          cwp = ciwpmcd(iplon,ig,lay) + clwpmcd(iplon,ig,lay)
		  ! (dmb 2012) the stop commands were removed because they aren't supported on the GPU
          if (cldfmc(iplon,ig,lay) .ge. cldmin .and. &
             (cwp .ge. cldmin .or. taucmc(iplon,ig,lay) .ge. cldmin)) then


            if(inflagd(iplon) .eq. 2) then
               radice = reicmcd(iplon, lay)

! Calculation of absorption coefficients due to ice clouds.
               if (ciwpmcd(iplon,ig,lay) .eq. 0.0 ) then
                  abscoice = 0.0 
				  
				 
               elseif (iceflagd(iplon) .eq. 0) then
                  
                  abscoice= absice0d(1) + absice0d(2)/radice

               elseif (iceflagd(iplon) .eq. 1) then
                    ncbands(iplon) = 5
                  ib = icb(ngb(ig))
                  abscoice = absice1d(1,ib) + absice1d(2,ib)/radice

! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns

               elseif (iceflagd(iplon) .eq. 2) then
                     ncbands(iplon) = 16
                     factor = (radice - 2. )/3. 
                     index = int(factor)
                     if (index .eq. 43) index = 42
                     fint = factor - float(index)
                     ib = ngb(ig)
                     abscoice = &
                         absice2d(index,ib) + fint * &
                         (absice2d(index+1,ib) - (absice2d(index,ib))) 
               
! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

               elseif (iceflagd(iplon) .eq. 3) then
                      ncbands(iplon) = 16
                     factor = (radice - 2. )/3. 
                     index = int(factor)
                     if (index .eq. 46) index = 45
                     fint = factor - float(index)
                     ib = ngb(ig)
                     abscoice= &
                         absice3d(index,ib) + fint * &
                         (absice3d(index+1,ib) - (absice3d(index,ib)))

! For iceflag=4 option, ice particle effective diameter is limited to 1.0 to 200.0 microns

               elseif (iceflagd(iplon) .eq. 4) then
                      ncbands(iplon) = 16
                     factor = radice
                     index = int(factor)
                     fint = factor - float(index)
                     ib = ngb(ig)
                     abscoice= &
                         absice4d(index,ib) + fint * &
                         (absice4d(index+1,ib) - (absice4d(index,ib)))
   
               endif
                  
! Calculation of absorption coefficients due to water clouds.
              if (liqflagd(iplon) .eq. 1) then
                  radliq = relqmcd(iplon, lay)
                  index = int(radliq - 1.5 )
                  if (index .eq. 0) index = 1
                  if (index .eq. 58) index = 57
                  fint = radliq - 1.5  - float(index)
                  ib = ngb(ig)
                  abscoliq = &
                        absliq1d(index,ib) + fint * &
                        (absliq1d(index+1,ib) - (absliq1d(index,ib)))
               endif

               taucmc(iplon,ig,lay) = ciwpmcd(iplon,ig,lay) * abscoice + &
                                clwpmcd(iplon,ig,lay) * abscoliq

            endif
         endif
        

#ifdef _CUDA
      endif
#else
      end do
      end do
      end do
#endif  


      end subroutine cldprmcg

      
      ! (dmb 2012) This subroutine allocates the module level arrays on the GPU
      subroutine allocateGPUcldprmcg(ncol, nlay, ngptlw)

	  integer , intent(in) :: nlay, ngptlw, ncol
      allocate( inflagd(ncol), iceflagd(ncol), liqflagd(ncol))
	  allocate( ciwpmcd(ncol,ngptlw,nlay+1))
	
	  allocate( relqmcd(ncol, nlay+1), reicmcd(ncol, nlay+1))
	  allocate( clwpmcd(ncol, ngptlw, nlay+1))
	
      
      
      end subroutine

      ! (dmb 2012) This subroutine deallocates any GPU arrays.
      subroutine deallocateGPUcldprmcg()

      deallocate( inflagd, iceflagd, liqflagd)
	  deallocate( ciwpmcd)
	
	  deallocate( relqmcd, reicmcd)
	  deallocate( clwpmcd)
	 
      
      
      end subroutine

      ! (dmb 2012) This subroutine copies input data from the CPU over to the GPU
      ! for use in the cldprmcg subroutine.
      subroutine copyGPUcldprmcg(inflag, iceflag, liqflag,&
                                 absice0, absice1, absice2, absice3, absice4, absliq1)
                                
      integer :: inflag(:), iceflag(:), liqflag(:)
	
     
       real , dimension(:) :: absice0
      real , dimension(:,:) :: absice1
      real , dimension(:,:) :: absice2
      real , dimension(:,:) :: absice3
      real , dimension(:,:) :: absice4
      real , dimension(:,:) :: absliq1
      
      inflagd = inflag
	  iceflagd = iceflag
	  liqflagd = liqflag


	  absice0d = absice0
	  absice1d = absice1
	  absice2d = absice2
	  absice3d = absice3
	  absice4d = absice4
	  absliq1d = absliq1
      
      
      
      end subroutine 


      end module gpu_rrtmg_lw_cldprmc
