   
   
#include "_gpudef.inc"
   
   module rrtmg_lw_rad


#ifdef _CUDA
      use cudafor

#endif

      use gpu_mcica_subcol_gen_lw


      use rrlw_vsn
      use gpu_rrtmg_lw_rtrnmc
      use gpu_rrtmg_lw_setcoef
      use gpu_rrtmg_lw_cldprmc


      use gpu_rrtmg_lw_taumol, only: taumolg, copyGPUTaumol
      use rrlw_cld, only: abscld1, absliq0, absliq1, &
            absice0, absice1, absice2, absice3, absice4
      use rrlw_wvn, only: ngb, ngs
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi, grav, avogad


      implicit none

      integer   _gpudev, allocatable :: ngbd(:)
      integer, allocatable _gpudev :: ncbandsd(:)
      integer,  allocatable _gpudev :: icbd(:)
      integer  , allocatable _gpudev :: icldlyr(:,:)
      real  _gpudev, allocatable :: fracsd(:,:,:)
      real  _gpudev, allocatable :: taug(:,:,:)

      real :: timings(10)

      !------------------------------------------------------------------
   contains
      !------------------------------------------------------------------
      ! (dmb 2012) Added the GPUFlag parameter
      ! GPUFlag = 0 to run on CPU, 1 to run on the GPU
      subroutine rrtmg_lw( &
            ncol ,nlay    ,icld    ,idrv    , &
            play    ,plev    ,tlay    ,tlev    ,tsfc    , & 
            h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr , &
            cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
            inflglw ,iceflglw,liqflglw,cldfrac , &
            tauc ,ciwp ,clwp ,rei ,rel , &
            tauaer  , zm, cloudMH, cloudHH, &
            uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, &
            duflx_dt, duflxc_dt, cloudFlag, &
            olrb06, dolrb06_dt, olrb09, dolrb09_dt, &
            olrb10, dolrb10_dt, olrb11, dolrb11_dt, &
            dyofyr, alat, numCPUs, partition_size)
         ! -------- Description --------

         ! This program is the driver subroutine for RRTMG_LW, the AER LW radiation 
         ! model for application to GCMs, that has been adapted from RRTM_LW for
         ! improved efficiency.
         !
         ! NOTE: The call to RRTMG_LW_INI should be moved to the GCM initialization
         !  area, since this has to be called only once. 
         !
         ! This routine:
         !    a) calls INATM to read in the atmospheric profile from GCM;
         !       all layering in RRTMG is ordered from surface to toa. 
         !    b) calls CLDPRMC to set cloud optical depth for McICA based 
         !       on input cloud properties 
         !    c) calls SETCOEF to calculate various quantities needed for 
         !       the radiative transfer algorithm
         !    d) calls TAUMOL to calculate gaseous optical depths for each 
         !       of the 16 spectral bands
         !    e) calls RTRNMC (for both clear and cloudy profiles) to perform the
         !       radiative transfer calculation using McICA, the Monte-Carlo 
         !       Independent Column Approximation, to represent sub-grid scale 
         !       cloud variability
         !    f) passes the necessary fluxes and cooling rates back to GCM
         !
         ! Two modes of operation are possible:
         !     The mode is chosen by using either rrtmg_lw.nomcica.f90 (to not use
         !     McICA) or rrtmg_lw.f90 (to use McICA) to interface with a GCM. 
         !
         !    1) Standard, single forward model calculation (imca = 0)
         !    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
         !       JC, 2003) method is applied to the forward model calculation (imca = 1)
         !
         ! This call to RRTMG_LW must be preceeded by a call to the module
         !     mcica_subcol_gen_lw.f90 to run the McICA sub-column cloud generator,
         !     which will provide the cloud physical or cloud optical properties
         !     on the RRTMG quadrature point (ngpt) dimension.
         !     Two random number generators are available for use when imca = 1.
         !     This is chosen by setting flag irnd on input to mcica_subcol_gen_lw.
         !     1) KISSVEC (irnd = 0)
         !     2) Mersenne-Twister (irnd = 1)
         !
         ! Two methods of cloud property input are possible:
         !     Cloud properties can be input in one of two ways (controlled by input 
         !     flags inflglw, iceflglw, and liqflglw; see text file rrtmg_lw_instructions
         !     and subroutine rrtmg_lw_cldprmc.f90 for further details):
         !
         !    1) Input cloud fraction and cloud optical depth directly (inflglw = 0)
         !    2) Input cloud fraction and cloud physical properties (inflglw = 1 or 2);  
         !       cloud optical properties are calculated by cldprmc or cldprmc based
         !       on input settings of iceflglw and liqflglw.  Ice particle size provided
         !       must be appropriately defined for the ice parameterization selected. 
         !
         ! One method of aerosol property input is possible:
         !     Aerosol properties can be input in only one way (controlled by input 
         !     flag iaer; see text file rrtmg_lw_instructions for further details):
         !
         !    1) Input aerosol optical depth directly by layer and spectral band (iaer=10);
         !       band average optical depth at the mid-point of each spectral band.
         !       RRTMG_LW currently treats only aerosol absorption;
         !       scattering capability is not presently available.
         !
         ! The optional calculation of the change in upward flux as a function of surface 
         ! temperature is available (controlled by input flag idrv).  This can be utilized 
         ! to approximate adjustments to the upward flux profile caused only by a change in 
         ! surface temperature between full radiation calls.  This feature uses the pre-
         ! calculated derivative of the Planck function with respect to surface temperature. 
         !
         !    1) Normal forward calculation for the input profile (idrv=0)
         !    2) Normal forward calculation with optional calculation of the change
         !       in upward flux as a function of surface temperature for clear sky
         !       and total sky flux.  Flux partial derivatives are provided in arrays
         !       duflx_dt and duflxc_dt for total and clear sky.  (idrv=1)
         !
         !
         ! ------- Modifications -------
         !
         ! This version of RRTMG_LW has been modified from RRTM_LW to use a reduced 
         ! set of g-points for application to GCMs.  
         !
         !-- Original version (derived from RRTM_LW), reduction of g-points, other
         !   revisions for use with GCMs.  
         !     1999: M. J. Iacono, AER, Inc.
         !-- Adapted for use with NCAR/CAM.
         !     May 2004: M. J. Iacono, AER, Inc.
         !-- Revised to add McICA capability. 
         !     Nov 2005: M. J. Iacono, AER, Inc.
         !-- Conversion to F90 formatting for consistency with rrtmg_sw.
         !     Feb 2007: M. J. Iacono, AER, Inc.
         !-- Modifications to formatting to use assumed-shape arrays.
         !     Aug 2007: M. J. Iacono, AER, Inc.
         !-- Modified to add longwave aerosol absorption.
         !     Apr 2008: M. J. Iacono, AER, Inc.
         !-- Added capability to calculate derivative of upward flux wrt surface temperature. 
         !     Nov 2009: M. J. Iacono, E. J. Mlawer, AER, Inc.
         !-- Added capability to run on GPU
         !     Aug 2012: David Berthiaume, AER, Inc.
         !-- Added alat and dyofyr to echo cloud-overlap scheme in RRTMGPU_SW
         !     Dec 2013: Matt Thompson, NASA/GMAO
         !-- Added numCPUs to allow for multple CPUs per GPU.
         !     Jan 2014: Matt Thompson, NASA/GMAO
         !-- Added new routine rtrnzero to help speed up zeroing on CPUs.
         !     Mar 2014: Matt Thompson, NASA/GMAO
         !-- Added band 10 OLR and d/dT for water vapor channel simulation
         !     Mar 2019: Peter Norris, NASA/GMAO
         ! --------- Modules ----------

         use parrrtm, only : nbndlw, ngptlw, maxxsec, mxmol, nbndlw
         use rrlw_con, only: fluxfac, heatfac, oneminus, pi
         use rrlw_wvn, only: ng, ngb, nspa, nspb, wavenum1, wavenum2, delwave

         ! ------- Declarations -------

         ! integer , parameter:: maxlay = 203
         ! integer , parameter:: mxmol = 38


         ! ----- Input -----
         ! Note: All volume mixing ratios are in dimensionless units of mole fraction obtained
         ! by scaling mass mixing ratio (g/g) with the appropriate molecular weights (g/mol) 
         integer , intent(in) :: ncol             ! Number of horizontal columns
         integer , intent(in) :: nlay             ! Number of model layers
         integer , intent(inout) :: icld          ! Cloud overlap method
         !    0: Clear only
         !    1: Random
         !    2: Maximum/random
         !    3: Maximum
         integer , intent(in) :: idrv             ! Flag for calculation of dFdT, the change
         !    in upward flux as a function of 
         !    surface temperature [0=off, 1=on]
         !    0: Normal forward calculation
         !    1: Normal forward calculation with
         !       duflx_dt and duflxc_dt output

         integer , intent(in) :: cloudMH, cloudHH ! cloud layer heights for cloudFlag
         real ,    intent(in) :: play(:,:)        ! Layer pressures (hPa, mb)
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: plev(:,0:)       ! Interface pressures (hPa, mb)
         !    Dimensions: (ncol,nlay+1)
         real ,    intent(in) :: tlay(:,:)        ! Layer temperatures (K)
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: tlev(:,0:)       ! Interface temperatures (K)
         !    Dimensions: (ncol,nlay+1)
         real ,    intent(in) :: tsfc(:)          ! Surface temperature (K)
         !    Dimensions: (ncol)
         real ,    intent(in) :: h2ovmr(:,:)      ! H2O volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: o3vmr(:,:)       ! O3 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: co2vmr(:,:)      ! CO2 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: ch4vmr(:,:)      ! Methane volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: n2ovmr(:,:)      ! Nitrous oxide volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: o2vmr(:,:)       ! Oxygen volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: cfc11vmr(:, :)   ! CFC11 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: cfc12vmr(:, :)   ! CFC12 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: cfc22vmr(:, :)   ! CFC22 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: ccl4vmr(:, :)    ! CCL4 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: emis(:, :)       ! Surface emissivity
         !    Dimensions: (ncol,nbndlw)

         integer , intent(in) :: inflglw          ! Flag for cloud optical properties
         integer , intent(in) :: iceflglw         ! Flag for ice particle specification
         integer , intent(in) :: liqflglw         ! Flag for liquid droplet specification

         real ,    intent(in) :: cldfrac(:,:)     ! Cloud fraction
         !    Dimensions: (ngptlw,ncol,nlay)
         real ,    intent(in) :: ciwp(:,:)        ! In-cloud ice water path (g/m2)
         !    Dimensions: (ngptlw,ncol,nlay)
         real ,    intent(in) :: clwp(:,:)        ! In-cloud liquid water path (g/m2)
         !    Dimensions: (ngptlw,ncol,nlay)
         real ,    intent(in) :: rei(:,:)         ! Cloud ice particle effective size (microns)
         !    Dimensions: (ncol,nlay)
         ! specific definition of reicmcl depends on setting of iceflglw:
         ! iceflglw = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
         !               r_ec must be >= 10.0 microns
         ! iceflglw = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
         !               r_ec range is limited to 13.0 to 130.0 microns
         ! iceflglw = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
         !               r_k range is limited to 5.0 to 131.0 microns
         ! iceflglw = 3: generalized effective size, dge, (Fu, 1996),
         !               dge range is limited to 5.0 to 140.0 microns
         !               [dge = 1.0315 * r_ec]
         real ,    intent(in) :: rel(:, :)        ! Cloud water drop effective radius (microns)
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: tauc(:, :, :)    ! In-cloud optical depth
         !    Dimensions: (ngptlw,ncol,nlay)

         real ,    intent(in) :: tauaer(:,:,:)    ! aerosol optical depth
         !   at mid-point of LW spectral bands
         !    Dimensions: (ncol,nlay,nbndlw)

         real ,    intent(in) :: zm(:,:)          ! Heights of level midpoints
         !    Dimensions(ncol,nlay)

         real ,    intent(in) :: alat(:)          ! Latitude of column
         !    Dimensions(ncol)
         integer , intent(in) :: dyofyr           ! Day of the year

         integer , intent(in) :: numCPUs          ! Number of cores 
         integer , intent(in) :: partition_size   ! Partition size

         ! ----- Output -----

         real , intent(out) :: uflx(:,:)          ! Total sky longwave upward flux (W/m2)
         !    Dimensions: (ncol,nlay+1)
         real , intent(out) :: dflx(:,:)          ! Total sky longwave downward flux (W/m2)
         !    Dimensions: (ncol,nlay+1)
         real , intent(out) :: hr(:,:)            ! Total sky longwave radiative heating rate (K/d)
         !    Dimensions: (ncol,nlay)
         real , intent(out) :: uflxc(:,:)         ! Clear sky longwave upward flux (W/m2)
         !    Dimensions: (ncol,nlay+1)
         real , intent(out) :: dflxc(:,:)         ! Clear sky longwave downward flux (W/m2)
         !    Dimensions: (ncol,nlay+1)
         real , intent(out) :: hrc(:,:)           ! Clear sky longwave radiative heating rate (K/d)
         !    Dimensions: (ncol,nlay)

         ! ----- Optional Output -----
         real , intent(out), optional :: duflx_dt(:,:)     
         ! change in upward longwave flux (w/m2/K)
         ! with respect to surface temperature
         !    Dimensions: (ncol,nlay)
         real , intent(out), optional :: duflxc_dt(:,:)    
         ! change in clear sky upward longwave flux (w/m2/K)
         ! with respect to surface temperature
         !    Dimensions: (ncol,nlay)
         integer , intent(out), optional :: cloudFlag(:,:)

         real, intent(out), dimension(:), optional :: olrb06, dolrb06_dt
         real, intent(out), dimension(:), optional :: olrb09, dolrb09_dt
         real, intent(out), dimension(:), optional :: olrb10, dolrb10_dt
         real, intent(out), dimension(:), optional :: olrb11, dolrb11_dt
         ! OLR for bands 9-11 and temperature derivatives (W/m2, W/m2/K)
         !    Dimensions: (ncol)

         integer  :: pncol
         integer  :: colstart
         integer  :: cn, ns, i, np, mns
         real :: minmem
         integer :: hetflag
         integer :: numDevices, err, numCPUsPerGPU

         integer :: numThreads
         ! Cuda device information
#ifdef _CUDA
         type(cudadeviceprop) :: prop
#endif
         ! store the available device global and constant memory
         real gmem, cmem
         real t1,t2

#ifdef _CUDA

         err = cudaGetDeviceProperties( prop, 0)
         gmem = prop%totalGlobalMem
         !print *, "total GPU global memory is ", gmem / (1024.0*1024.0) , "MB"

         err = cudaGetDeviceCount(numDevices)
         !print *, "total number of GPUs is ", numDevices
         !print *, "total number of CPUs is ", numCPUs

         numCPUsPerGPU = ceiling( real(numCPUs) / real(numDevices) )
         !print *, "number of CPUs per GPU is", numCPUsPerGPU

#endif

         ! (dmb 2012) Here we calculate the number of groups to partition
         ! the inputs.

         ! determine the minimum GPU memory
         ! force the GPUFlag off if there are no devices available

#ifdef _CUDA
         minmem = gmem/(real(numCPUsPerGPU))

         ! MAT At high core count, this seems to not be sufficient and 
         !     it can cause an issue with too little memory per CPU on 
         !     the GPU.  So, we reduce minmem, to cause cn below to get 
         !     smaller. The factor is...a guess, really.

         minmem = minmem / 3.0

         ! use the available memory to determine the minumum number 
         ! of steps that will be required.
         ! We use 1000 profiles per available GB as a conservative 
         ! lower bound.
         cn = minmem * 1000 / (1024**3)
#else
         ! MAT For reasons as yet unknown, partitioning on CPU leads to 
         !     different answers than when you do not. For now, set this
         !     value absurdly high so we never partition. 33000 is chosen
         !     as in case someone wants to run C180 on 1x6, this is
         !     greater than 180*180=32400
         minmem = 33.0 * (1024.0**3)

         if (partition_size < 0) then

            ! use the available memory to determine the minumum number 
            ! of steps that will be required.
            ! We use 1000 profiles per available GB as a conservative 
            ! lower bound.
            cn = minmem * 1000 / (1024**3)
            !print *, 'cn with calc: ', cn
            ! with device emulation (for debugging) make sure there is a lower
            ! limit to the number of supported columns
            !if (cn < 500) then 
            !cn = 500 
            !end if
         else
            cn = partition_size
         end if
             
         ! set the number of 'devices' to the available number of CPUs
#endif
         !print *, "available working memory is ", int(minmem / (1024*1024)) , " MB"

         !print *, "number of profiles per partition is ", cn
         ns = ceiling( real(ncol) / real(cn) )

         !print *, "number of steps is set to ", ns

         call cpu_time(t1)

         do  i = 1, ns 
            !print *, "running partition ", i , " of ", ns
            !print *, "running partition ", i , " of ", ns, " ncol ", ncol , " colstart ", (i-1)*cn + 1, " pncol ", min(cn, ncol - (i-1)*cn)


            call rrtmg_lw_part( ns, ncol, (i-1)*cn + 1, min(cn, ncol - (i-1)*cn), nlay, icld, idrv,&
                  play    ,plev    ,tlay    ,tlev    ,tsfc    , & 
                  h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr , &
                  cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
                  inflglw ,iceflglw,liqflglw,cldfrac , &
                  tauc ,ciwp ,clwp ,rei ,rel , &
                  tauaer  , zm, cloudMH, cloudHH, &
                  uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, &
                  duflx_dt, duflxc_dt, cloudFlag, &
                  olrb06, dolrb06_dt, olrb09, dolrb09_dt, &
                  olrb10, dolrb10_dt, olrb11, dolrb11_dt, &
                  dyofyr,alat)    
         end do


         call cpu_time(t2)
         !print *, "------------------------------------------------"
         !print *, "TOTAL RUN TIME IS   ", t2-t1

      end subroutine rrtmg_lw




      subroutine rrtmg_lw_part &
            (npart, ncol , colstart, pncol ,nlay    ,icld    ,idrv    , &
            play    ,plev    ,tlay    ,tlev    ,tsfc    , & 
            h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr , &
            cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
            inflglw ,iceflglw,liqflglw,cldfrac , &
            tauc ,ciwp ,clwp ,rei ,rel , &
            tauaer  , zm, cloudMH, cloudHH, &
            uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, &
            duflx_dt, duflxc_dt, cloudFlag, &
            olrb06, dolrb06_dt, olrb09, dolrb09_dt, &
            olrb10, dolrb10_dt, olrb11, dolrb11_dt, &
            dyofyr,alat)


         use gpu_mcica_subcol_gen_lw, only: mcica_subcol_lwg, mcica_cldc, generate_stochastic_cloudsg

         use parrrtm, only : nbndlw, ngptlw, maxxsec, mxmol, nbndlw, nmol
         use rrlw_con, only: fluxfac, heatfac, oneminus, pi
         use rrlw_wvn, only: ng, ngb, nspa, nspb, wavenum1, wavenum2, delwave, ixindx



         ! ----- Input -----
         ! Note: All volume mixing ratios are in dimensionless units of mole fraction obtained
         ! by scaling mass mixing ratio (g/g) with the appropriate molecular weights (g/mol) 
         integer , intent(in) :: npart
         integer , intent(in) :: ncol              ! Number of horizontal columns
         integer , intent(in) :: nlay              ! Number of model layers
         integer , intent(inout) :: icld           ! Cloud overlap method
         !    0: Clear only
         !    1: Random
         !    2: Maximum/random
         !    3: Maximum
         integer , intent(in) :: idrv              ! Flag for calculation of dFdT, the change
         !    in upward flux as a function of 
         !    surface temperature [0=off, 1=on]
         !    0: Normal forward calculation
         !    1: Normal forward calculation with
         !       duflx_dt and duflxc_dt output

         real ,    intent(in) :: play(:,:)        ! Layer pressures (hPa, mb)
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: plev(:,0:)       ! Interface pressures (hPa, mb)
         !    Dimensions: (ncol,nlay+1)
         real ,    intent(in) :: tlay(:,:)        ! Layer temperatures (K)
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: tlev(:,0:)       ! Interface temperatures (K)
         !    Dimensions: (ncol,nlay+1)
         real ,    intent(in) :: tsfc(:)          ! Surface temperature (K)
         !    Dimensions: (ncol)
         real ,    intent(in) :: h2ovmr(:,:)      ! H2O volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: o3vmr(:,:)       ! O3 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: co2vmr(:,:)      ! CO2 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: ch4vmr(:,:)      ! Methane volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: n2ovmr(:,:)      ! Nitrous oxide volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: o2vmr(:,:)       ! Oxygen volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: cfc11vmr(:, :)   ! CFC11 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: cfc12vmr(:, :)   ! CFC12 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: cfc22vmr(:, :)   ! CFC22 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: ccl4vmr(:, :)    ! CCL4 volume mixing ratio
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: emis(:, :)       ! Surface emissivity
         !    Dimensions: (ncol,nbndlw)

         integer , intent(in) :: inflglw          ! Flag for cloud optical properties
         integer , intent(in) :: iceflglw         ! Flag for ice particle specification
         integer , intent(in) :: liqflglw         ! Flag for liquid droplet specification

         integer , intent(in) :: cloudMH, cloudHH ! cloud layer heights for cloudFlag
         real ,    intent(in) :: cldfrac(:,:)     ! Cloud fraction
         !    Dimensions: (ngptlw,ncol,nlay)
         real ,    intent(in) :: ciwp(:,:)        ! In-cloud ice water path (g/m2)
         !    Dimensions: (ngptlw,ncol,nlay)
         real ,    intent(in) :: clwp(:,:)        ! In-cloud liquid water path (g/m2)
         !    Dimensions: (ngptlw,ncol,nlay)
         real ,    intent(in) :: rei(:,:)         ! Cloud ice particle effective size (microns)
         !    Dimensions: (ncol,nlay)
         ! specific definition of reicmcl depends on setting of iceflglw:
         ! iceflglw = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
         !               r_ec must be >= 10.0 microns
         ! iceflglw = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
         !               r_ec range is limited to 13.0 to 130.0 microns
         ! iceflglw = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
         !               r_k range is limited to 5.0 to 131.0 microns
         ! iceflglw = 3: generalized effective size, dge, (Fu, 1996),
         !               dge range is limited to 5.0 to 140.0 microns
         !               [dge = 1.0315 * r_ec]
         real ,    intent(in) :: rel(:, :)        ! Cloud water drop effective radius (microns)
         !    Dimensions: (ncol,nlay)
         real ,    intent(in) :: tauc(:, :,:)     ! In-cloud optical depth
         !    Dimensions: (ngptlw,ncol,nlay)

         real ,    intent(in) :: tauaer(:,:,:)    ! aerosol optical depth
         !   at mid-point of LW spectral bands
         !    Dimensions: (ncol,nlay,nbndlw)
         real ,    intent(in) :: zm(:, :)         ! Heights of level midpoints
         !    Dimensions(ncol,nlay)

         real ,    intent(in) :: alat(:)          ! Latitude of column
         !    Dimensions(ncol)
         integer , intent(in) :: dyofyr           ! Day of the year


         integer , intent(in) :: pncol
         integer , intent(in) :: colstart

         ! ----- Output -----

         real , intent(out) :: uflx(:,:)          ! Total sky longwave upward flux (W/m2)
         !    Dimensions: (ncol,nlay+1)
         real , intent(out) :: dflx(:,:)          ! Total sky longwave downward flux (W/m2)
         !    Dimensions: (ncol,nlay+1)
         real , intent(out) :: hr(:,:)            ! Total sky longwave radiative heating rate (K/d)
         !    Dimensions: (ncol,nlay)
         real , intent(out) :: uflxc(:,:)         ! Clear sky longwave upward flux (W/m2)
         !    Dimensions: (ncol,nlay+1)
         real , intent(out) :: dflxc(:,:)         ! Clear sky longwave downward flux (W/m2)
         !    Dimensions: (ncol,nlay+1)
         real , intent(out) :: hrc(:,:)           ! Clear sky longwave radiative heating rate (K/d)
         !    Dimensions: (ncol,nlay)

         ! ----- Optional Output -----
         real , intent(out), optional :: duflx_dt(:,:)     
         ! change in upward longwave flux (w/m2/K)
         ! with respect to surface temperature
         !    Dimensions: (ncol,nlay)
         real , intent(out), optional :: duflxc_dt(:,:)    
         ! change in clear sky upward longwave flux (w/m2/K)
         ! with respect to surface temperature
         !    Dimensions: (ncol,nlay)
         integer , intent(out), optional :: cloudFlag(:,:)

         real, intent(out), dimension(:), optional :: olrb06, dolrb06_dt
         real, intent(out), dimension(:), optional :: olrb09, dolrb09_dt
         real, intent(out), dimension(:), optional :: olrb10, dolrb10_dt
         real, intent(out), dimension(:), optional :: olrb11, dolrb11_dt
         ! OLR for bands 6 & 9-11 and temperature derivatives (W/m2, W/m2/K)
         !    Dimensions: (ncol)

         real  _gpudeva :: cldfmcd(:,:,:)         ! layer cloud fraction [mcica]
         !    Dimensions: (ngptlw,nlayers)


!GPU These are needed to allocate space for temporary automatic local
!GPU arrays on the GPU. On the CPU this ifndef converts them to nlay

#ifndef GPU_MAXLEVS
#define GPU_MAXLEVS nlay
#endif


         ! ----- Local -----

         ! Control
         integer(kind=4) :: nlayers            ! total number of layers
         integer(kind=4) :: istart             ! beginning band of calculation
         integer(kind=4) :: iend               ! ending band of calculation
         integer(kind=4) :: iout               ! output option flag (inactive)
         integer  :: iaer                      ! aerosol option flag
         integer(kind=4) :: iplon              ! column loop index
         integer  :: imca                      ! flag for mcica [0=off, 1=on]
         integer  :: ims                       ! value for changing mcica permute seed
         integer  :: k                         ! layer loop index
         integer  :: ig                        ! g-point loop index
         real  :: t1, t2
         ! Atmosphere
         real  :: pavel(pncol,GPU_MAXLEVS)          ! layer pressures (mb) 
         real  :: tavel(pncol,GPU_MAXLEVS)          ! layer temperatures (K)
         real  :: pz(pncol,0:GPU_MAXLEVS)           ! level (interface) pressures (hPa, mb)
         real  :: tz(pncol,0:GPU_MAXLEVS)           ! level (interface) temperatures (K)
         real  :: tbound(pncol)                ! surface temperature (K)
         real  :: coldry(pncol,GPU_MAXLEVS)         ! dry air column density (mol/cm2)
         real  :: wbrodl(pncol,GPU_MAXLEVS)         ! broadening gas column density (mol/cm2)
         !real  :: wkl(pncol,mxmol,GPU_MAXLEVS+1)      ! molecular amounts (mol/cm-2)
         !real  :: wx(pncol,maxxsec,GPU_MAXLEVS+1)     ! cross-section amounts (mol/cm-2)
         real  :: pwvcm(pncol)                 ! precipitable water vapor (cm)
         real  :: semiss(pncol,nbndlw)         ! lw surface emissivity


         ! Atmosphere/clouds - cldprop
         integer  :: ncbands(pncol)                   ! number of cloud spectral bands
         integer  :: inflag(pncol)              ! flag for cloud property method
         integer  :: iceflag(pncol)             ! flag for ice cloud properties
         integer  :: liqflag(pncol)             ! flag for liquid cloud properties



         ! Output
         !real  :: totuflux(pncol,0:GPU_MAXLEVS)     ! upward longwave flux (w/m2)
         !real  :: totdflux(pncol,0:GPU_MAXLEVS)     ! downward longwave flux (w/m2)
         !real  :: fnet(pncol,0:GPU_MAXLEVS)         ! net longwave flux (w/m2)
         !real  :: htr(pncol,0:GPU_MAXLEVS)          ! longwave heating rate (k/day)
         !real  :: totuclfl(pncol,0:GPU_MAXLEVS)     ! clear sky upward longwave flux (w/m2)
         !real  :: totdclfl(pncol,0:GPU_MAXLEVS)     ! clear sky downward longwave flux (w/m2)
         !real  :: fnetc(pncol,0:GPU_MAXLEVS)        ! clear sky net longwave flux (w/m2)
         !real  :: htrc(pncol,0:GPU_MAXLEVS)         ! clear sky longwave heating rate (k/day)
         !real  :: dtotuflux_dt(pncol,0:GPU_MAXLEVS) ! change in upward longwave flux (w/m2/k)
         ! with respect to surface temperature
         !real  :: dtotuclfl_dt(pncol,0:GPU_MAXLEVS) ! change in clear sky upward longwave flux (w/m2/k)
         ! with respect to surface temperature
         !real  ::  curad(pncol,ngptlw,0:GPU_MAXLEVS)     ! upward longwave flux (w/m2)
         !real  ::   cdrad(pncol,ngptlw,0:GPU_MAXLEVS)     ! downward longwave flux (w/m2)
         !real  ::   cclrurad(pncol,ngptlw,0:GPU_MAXLEVS)     ! clear sky upward longwave flux (w/m2)
         !real  ::   cclrdrad(pncol,ngptlw,0:GPU_MAXLEVS)     ! clear sky downward longwave flux (w/m2)
         !real  :: olrb10(pncol)      ! TOA OLR in band10 (W/m2)
         !real  :: dolrb10_dt(pncol)  ! change in TOA OLR in band10 (W/m2/K)

         real  :: cldfracq(pncol,GPU_MAXLEVS)     ! Cloud fraction
         !    Dimensions: (ngptlw,ncol,nlay)
         real  :: ciwpq(pncol,GPU_MAXLEVS)     ! In-cloud ice water path (g/m2)
         !    Dimensions: (ngptlw,ncol,nlay)
         real  :: clwpq(pncol,GPU_MAXLEVS)     ! In-cloud liquid water path (g/m2)
         !    Dimensions: (ngptlw,ncol,nlay)
         real  :: reiq(pncol,GPU_MAXLEVS)       ! Cloud ice particle effective size (microns)
         !    Dimensions: (ncol,nlay)
         ! specific definition of reicmcl depends on setting of iceflglw:
         ! iceflglw = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
         !               r_ec must be >= 10.0 microns
         ! iceflglw = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
         !               r_ec range is limited to 13.0 to 130.0 microns
         ! iceflglw = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
         !               r_k range is limited to 5.0 to 131.0 microns
         ! iceflglw = 3: generalized effective size, dge, (Fu, 1996),
         !               dge range is limited to 5.0 to 140.0 microns
         !               [dge = 1.0315 * r_ec]
         real  :: relq(pncol, GPU_MAXLEVS)       ! Cloud water drop effective radius (microns)
         !    Dimensions: (ncol,nlay)
         real  :: taucq(pncol, nbndlw, GPU_MAXLEVS)     ! In-cloud optical depth
         !    Dimensions: (ngptlw,ncol,nlay)



         real  :: zmp(pncol, GPU_MAXLEVS)

         real  :: alatp(pncol)

         integer  :: icb(16)

         ! local looping variables
         integer :: i,j,kk, piplon

         ! cuda return code
         integer :: ierr
         ! cuda thread and grid block dimensions
#ifdef _CUDA
         type(dim3) :: dimGrid, dimBlock
#endif

         real , dimension(16) :: a0 =(/ 1.66 ,  1.55 ,  1.58 ,  1.66 , &
               1.54 , 1.454 ,  1.89 ,  1.33 , &
               1.668 ,  1.66 ,  1.66 ,  1.66 , &
               1.66 ,  1.66 ,  1.66 ,  1.66  /)
         real , dimension(16) :: a1=(/ 0.00 ,  0.25 ,  0.22 ,  0.00 , &
               0.13 , 0.446 , -0.10 ,  0.40 , &
               -0.006 ,  0.00 ,  0.00 ,  0.00 , &
               0.00 ,  0.00 ,  0.00 ,  0.00  /)
         real , dimension(16) :: a2 =(/ 0.00 , -12.0 , -11.7 ,  0.00 , &
               -0.72 ,-0.243 ,  0.19 ,-0.062 , &
               0.414 ,  0.00 ,  0.00 ,  0.00 , &
               0.00 ,  0.00 ,  0.00 ,  0.00  /)
         real , parameter :: amd = 28.9660     ! Effective molecular weight of dry air (g/mol)
         real , parameter :: amw = 18.0160     ! Molecular weight of water vapor (g/mol)

         ! (dmb 2012) these arrays were moved to the main routine so that we can bypass some of the 
         ! inatm inefficiencies when running on the GPU
         real , parameter :: amdw = 1.607793   ! Molecular weight of dry air / water vapor
         real , parameter :: amdc = 0.658114   ! Molecular weight of dry air / carbon dioxide
         real , parameter :: amdo = 0.603428   ! Molecular weight of dry air / ozone
         real , parameter :: amdm = 1.805423   ! Molecular weight of dry air / methane
         real , parameter :: amdn = 0.658090   ! Molecular weight of dry air / nitrous oxide
         real , parameter :: amdo2 = 0.905140  ! Molecular weight of dry air / oxygen
         real , parameter :: amdc1 = 0.210852  ! Molecular weight of dry air / CFC11
         real , parameter :: amdc2 = 0.239546  ! Molecular weight of dry air / CFC12
         real  _gpudeva ::  zmd(:,:)
         real  :: amm, amttl, wvttl, wvsh, summol  
         integer  :: isp, l, ix, n, imol, ib   ! Loop indices
         integer, save :: counter =0
         real  :: btemp
         integer  :: cloudFlagq(pncol, 4)
         integer _gpudev :: pncold, nlayd, icldd, dyofyrd
         real  _gpudeva ::  alatd(:)


         !
         ! Initializations
         icb(:) = (/  1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5 /)

         oneminus = 1.  - 1.e-6 
         pi = 2.  * asin(1. )
         fluxfac = pi * 2.e4                   ! orig:   fluxfac = pi * 2.d4  
         istart = 1
         iend = 16
         iout = 0
         ims = 1
         pncold = pncol
         nlayd = nlay



         cldfracq(:,1:nlay) = cldfrac(colstart:(colstart+pncol-1), 1:nlay)
         ciwpq(:,1:nlay) = ciwp(colstart:(colstart+pncol-1), 1:nlay)
         clwpq(:,1:nlay) = clwp(colstart:(colstart+pncol-1), 1:nlay)
         reiq(:,1:nlay) = rei(colstart:(colstart+pncol-1), 1:nlay)
         relq(:,1:nlay) = rel(colstart:(colstart+pncol-1), 1:nlay)
         taucq(:,:,1:nlay) = tauc(colstart:(colstart+pncol-1), :, 1:nlay)
         zmp(:,1:nlay) = zm(colstart:(colstart+pncol-1), 1:nlay)


         allocate( cldfmcd(pncol, ngptlw, nlay+1))
         allocate( ngbd(140) )


         allocate( icbd(16))
         allocate( ncbandsd(pncol))
         allocate( icldlyr(pncol, nlay+1))



         call allocateGPUcldprmcg(pncol, nlay, ngptlw)
         call allocateGPUrtrnmcg(pncol, nlay, ngptlw, idrv)

         ngbd = ngb
         ngsd = ngs
         icldd = icld
         dyofyrd = dyofyr




         ! Set imca to select calculation type:
         !  imca = 0, use standard forward model calculation
         !  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability

         ! *** This version uses McICA (imca = 1) ***

         ! Set icld to select of clear or cloud calculation and cloud overlap method  
         ! icld = 0, clear only
         ! icld = 1, with clouds using random cloud overlap
         ! icld = 2, with clouds using maximum/random cloud overlap
         ! icld = 3, with clouds using maximum cloud overlap (McICA only)
         if (icld.lt.0.or.icld.gt.4) icld = 2

         ! Set iaer to select aerosol option
         ! iaer = 0, no aerosols
         ! icld = 10, input total aerosol optical depth (tauaer) directly
         iaer = 10

         ! Call model and data initialization, compute lookup tables, perform
         ! reduction of g-points from 256 to 140 for input absorption coefficient 
         ! data and other arrays.
         !
         ! In a GCM this call should be placed in the model initialization
         ! area, since this has to be called only once.  
         ! call rrtmg_lw_ini(cpdair)

         !     call rrtmg_lw_ini(1.004 )
         !  This is the main longitude/column loop within RRTMG.
         !  Prepare atmospheric profile from GCM for use in RRTMG, and define
         !  other input parameters.  



         ! (dmb 2012)

         nlayers = nlay

         call allocateGPUTaumol( pncol, nlayers, npart)

         allocate( fracsd( pncol, nlayers+1, ngptlw ))
         allocate( taug( pncol, nlayers+1, ngptlw ))
         allocate( zmd(pncol, nlay))
         allocate( alatd(pncol))
         tbound = tsfc(colstart:(colstart+pncol-1))
         pz(:,0:nlay) = plev(colstart:(colstart+pncol-1),0:nlay)
         tz(:,0:nlay) = tlev(colstart:(colstart+pncol-1),0:nlay)
         pavel(:,1:nlay) = play(colstart:(colstart+pncol-1),1:nlay)
         tavel(:,1:nlay) = tlay(colstart:(colstart+pncol-1),1:nlay)

         alatp = alat(colstart:(colstart+pncol-1))



         call copyGPUTaumolMol( colstart, pncol, nlayers, h2ovmr, co2vmr, o3vmr, n2ovmr, ch4vmr, &
               o2vmr, ccl4vmr, cfc11vmr, cfc12vmr, cfc22vmr, npart)



         zmd = zmp
         alatd = alatp

         call mcica_subcol_lwg(1, pncol, nlay, icld, counter, 0, pavel, cldfracq, ciwpq, &
               clwpq, taucq,ngbd, cldfmcd, ciwpmcd, clwpmcd, & 
               taucmcd,  cloudFlagq, cloudMH, cloudHH, zmd, alatd, dyofyr, 1)

         !  Generate the stochastic subcolumns of cloud optical properties for the longwave;
#ifdef _CUDA
         dimGrid = dim3( (ncol+255)/256,(140+1)/2, 1)
         dimBlock = dim3( 256,2,1)
#endif

         if (icld > 0) then
            call generate_stochastic_cloudsg _gpuchv (pncold, nlayd, icldd, &
                  ngbd, cldfmcd, clwpmcd, ciwpmcd, taucmcd, 1, zmd)
         end if

         do iplon = 1, pncol

            piplon = iplon + colstart - 1
            amttl = 0.0 
            wvttl = 0.0 
            do l = 1, nlayers

               amm = (1.  - h2ovmr(piplon,l)) * amd +h2ovmr(piplon,l) * amw            
               coldry(iplon, l) = (pz(iplon, l-1)-pz(iplon, l)) * 1.e3  * avogad / &
                     (1.e2  * grav * amm * (1.  + h2ovmr(piplon,l)))
            end do

            do l = 1, nlayers


               summol = co2vmr(piplon,l) + o3vmr(piplon,l) + n2ovmr(piplon,l) + ch4vmr(piplon,l) + o2vmr(piplon,l) 
               btemp = h2ovmr(piplon, l) * coldry(iplon, l)
               wbrodl(iplon, l) = coldry(iplon, l) * (1.  - summol)
               amttl = amttl + coldry(iplon, l)+btemp
               wvttl = wvttl + btemp
            enddo

            wvsh = (amw * wvttl) / (amd * amttl)
            pwvcm(iplon) = wvsh * (1.e3  * pz(iplon, 0)) / (1.e2  * grav)


            ! Transfer aerosol optical properties to RRTM variable;
            ! modify to reverse layer indexing here if necessary.

            if (icld .ge. 1) then 
               inflag(iplon) = inflglw
               iceflag(iplon) = iceflglw
               liqflag(iplon) = liqflglw

               ! Move incoming GCM cloud arrays to RRTMG cloud arrays.
               ! For GCM input, incoming reicmcl is defined based on selected ice parameterization (inflglw)

            endif
         enddo




         call mcica_cldc(pncol, icld, nlay, cloudmh, cloudhh, cloudflagq, cldfmcd)
         cloudFlag(colstart:(colstart+pncol-1), :) = cloudFlagq
         deallocate(zmd)
         deallocate(alatd)

         !  For cloudy atmosphere, use cldprmc to set cloud optical properties based on
         !  input cloud physical properties.  Select method based on choices described
         !  in cldprmc.  Cloud fraction, water path, liquid droplet and ice particle
         !  effective radius must be passed into cldprmc.  Cloud fraction and cloud
         !  optical depth are transferred to rrtmg_lw arrays in cldprmc.  

         ! If the GPU flag is active, then we call the GPU code.  Otherwise, call the CPU code  


         ! (dmb 2012) Copy the needed arrays over to the GPU for the cldprmc subroutine.
         call copyGPUcldprmcg( inflag, iceflag, liqflag,&
               absice0, absice1, absice2, absice3, absice4, absliq1 )


         ! copy common arrays over to the GPU
         icbd = icb
         a0d=a0
         a1d=a1
         a2d=a2
         delwaved=delwave
         relqmcd = relq
         reicmcd = reiq

         icldlyr = 0.0
         ! (dmb 2012) Allocate the arrays for the SetCoef and Taumol kernels
         call allocateGPUSetCoef( pncol, nlayers)

         ! (dmb 2012) Copy the needed data of to the GPU for the SetCoef and Taumol kernels


         call copyGPUTaumol( pavel, coldry, tauaer, pncol, colstart, nlay , npart)



         call copyGPUSetCoef( )
         ! (dmb 2012) Copy over additional common arrays 
         taveld = tavel
         tzd = tz
         tboundd = tbound
         wbroadd = wbrodl
         !wkld = wkl
         semissd = emis(colstart:(colstart+pncol-1),1:nbndlw)
         call copyToGPUref()
         call copyGPUrtrnmcg(pz, pwvcm, idrv)



         ! (dmb 2012) Here we configure the grids and blocks to run the cldpmcd kernel
         ! on the GPU.  I decided to keep the block dimensions to 16x16 to coincide with
         ! coalesced memory access when I am able to parition the profiles to multiples
         ! of 32.
#ifdef _CUDA
         dimGrid = dim3( (pncol+255)/256,(nlayers)/1, ngptlw)
         dimBlock = dim3( 256,1,1)
#endif
         !clwpmcd = 0
         !clwpmcd = clwpmc
         ! (dmb 2012) Call the cldprmcg kernel
         call cldprmcg _gpuchv (pncol, nlayers, cldfmcd, taucmcd,  ngbd, icbd, ncbandsd, icldlyr)

         ! synchronize the GPU with the CPU before taking timing results or passing data back to the CPU
#ifdef _CUDA   
         ierr = cudaThreadSynchronize() 
#endif



         ! Calculate information needed by the radiative transfer routine
         ! that is specific to this atmosphere, especially some of the 
         ! coefficients and indices needed to compute the optical depths
         ! by interpolating data from stored reference atmospheres. 




         ! (dmb 2012) Initialize the grid and block dimensions and call the setcoefg kernel
#ifdef _CUDA
         dimGrid = dim3( (pncol+255)/256,1, 1)
         dimBlock = dim3( 256,1,1)
#endif
         call setcoefg _gpuchv (pncol, nlayers, istart,idrv)
         ! (dmb 2012) end if GPU flag 




         !  Calculate the gaseous optical depths and Planck fractions for 
         !  each longwave spectral band.



         ! (dmb 2012) Call the taumolg subroutine.  This subroutine calls all of the individal taumol kernels.     
         call taumolg(1, pncol,nlayers, ngbd, taug, fracsd)





         ! Call the radiative transfer routine.
         ! Either routine can be called to do clear sky calculation.  If clouds
         ! are present, then select routine based on cloud overlap assumption
         ! to be used.  Clear sky calculation is done simultaneously.
         ! For McICA, RTRNMC is called for clear and cloudy calculations.

         ! MAT Added new routine rrtnzero to allow compiler to speed up 
         !     zeroing on CPUs.

#ifdef _CUDA
         ierr = cudaThreadSynchronize()
#endif 


#ifdef _CUDA    
         dimGrid = dim3( (pncol+255)/256, 70, 1)
         dimBlock = dim3( 256,2,1)
#endif    

         call rtrnzero _gpuchv (pncol,nlayers)

#ifdef _CUDA
         ierr = cudaThreadSynchronize()
#endif 


#ifdef _CUDA    
         dimGrid = dim3( (pncol+255)/256, 70, 1)
         dimBlock = dim3( 256,2,1)
#endif    

         call rtrnmcg _gpuchv (pncol,nlayers, istart, iend, iout, &
               ngbd, icldlyr, taug, fracsd, cldfmcd)

#ifdef _CUDA
         ierr = cudaThreadSynchronize() 
#endif




         ! sum up the results

         totufluxd = 0.0
         totdfluxd = 0.0
         totuclfld = 0.0
         totdclfld = 0.0
         htrd = 0.0
         htrcd = 0.0
         dtotuflux_dtd = 0.0
         dtotuclfl_dtd = 0.0
         olrb06d = 0.0
         dolrb06_dtd = 0.0
         olrb09d = 0.0
         dolrb09_dtd = 0.0
         olrb10d = 0.0
         dolrb10_dtd = 0.0
         olrb11d = 0.0
         dolrb11_dtd = 0.0

#ifdef CUDA
         dimGrid = dim3( (pncol+255)/256,nlayers+1,1)
         dimBlock = dim3( 256, 1, 1)
#endif
         ! (dmb 2012) Here we integrate across the g-point fluxes to arrive at total fluxes
         ! This functionality was factored out of the original rtrnmc routine so that I could
         ! parallelize across multiple dimensions.
         call rtrnadd _gpuchv (pncol, nlayers, ngptlw, idrv, ngbd)

#ifdef _CUDA
         ierr = cudaThreadSynchronize() 
         dimGrid = dim3( (pncol+255)/256,nlayers,1)
#endif


         ! (dmb 2012) Calculate the heating rates.
         call rtrnheatrates _gpuchv (pncol, nlayers)    

#ifdef CUDA
         ierr = cudaThreadSynchronize() 
#endif    
         ! copy the partition data back to the CPU
         uflx(colstart:(colstart+pncol-1), 1:(nlayers+1)) = totufluxd(:,0:nlayers)
         dflx(colstart:(colstart+pncol-1), 1:(nlayers+1)) = totdfluxd(:,0:nlayers)
         uflxc(colstart:(colstart+pncol-1), 1:(nlayers+1)) = totuclfld(:,0:nlayers)
         dflxc(colstart:(colstart+pncol-1), 1:(nlayers+1)) = totdclfld(:,0:nlayers)
         hr(colstart:(colstart+pncol-1), 1:(nlayers+1)) = htrd(:,0:nlayers)
         hrc(colstart:(colstart+pncol-1), 1:(nlayers+1)) = htrcd(:,0:nlayers)
         olrb06(colstart:(colstart+pncol-1)) = olrb06d(:)
         olrb09(colstart:(colstart+pncol-1)) = olrb09d(:)
         olrb10(colstart:(colstart+pncol-1)) = olrb10d(:)
         olrb11(colstart:(colstart+pncol-1)) = olrb11d(:)

         if (idrv .eq. 1) then

            duflx_dt(colstart:(colstart+pncol-1), 1:(nlayers+1)) = dtotuflux_dtd(:,0:nlayers)
            duflxc_dt(colstart:(colstart+pncol-1), 1:(nlayers+1)) = dtotuclfl_dtd(:,0:nlayers)
            dolrb06_dt(colstart:(colstart+pncol-1)) = dolrb06_dtd(:)
            dolrb09_dt(colstart:(colstart+pncol-1)) = dolrb09_dtd(:)
            dolrb10_dt(colstart:(colstart+pncol-1)) = dolrb10_dtd(:)
            dolrb11_dt(colstart:(colstart+pncol-1)) = dolrb11_dtd(:)

         end if






         !  Transfer up and down fluxes and heating rate to output arrays.
         !  Vertical indexing goes from bottom to top; reverse here for GCM if necessary.

         deallocate( cldfmcd)
         deallocate( icbd)
         deallocate( ncbandsd)
         deallocate( icldlyr)

         call deallocateGPUTaumol()
         deallocate( fracsd)
         deallocate( taug)
         deallocate( ngbd)
         call deallocateGPUcldprmcg()
         call deallocateGPUrtrnmcg(idrv)
         call deallocateGPUSetCoef( )

      end subroutine rrtmg_lw_part

   end module  rrtmg_lw_rad
