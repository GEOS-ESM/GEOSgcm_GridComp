!   $Id$

module soradmod

#ifdef _CUDA
   use cudafor
#else
   use sorad_constants, only: &
         wk_uv, zk_uv, ry_uv, &
         xk_ir, ry_ir,        &
         cah,   coa
   use rad_constants, only: &
         aig_uv, awg_uv, arg_uv, &
         aib_uv, awb_uv, arb_uv, &
         aib_nir, awb_nir, arb_nir, &
         aia_nir, awa_nir, ara_nir, &
         aig_nir, awg_nir, arg_nir, &
         caib, caif
   use gettau
#endif

   use MAPL_ConstantsMod, only: MAPL_R4, MAPL_R8, MAPL_GRAV

   implicit none

#ifndef _CUDA
   private

   public sorad
#endif

#ifdef _CUDA   

   ! Device inputs
   ! -------------

   real, allocatable, device :: cosz_dev(:)
   real, allocatable, device :: pl_dev(:,:)
   real, allocatable, device :: ta_dev(:,:)
   real, allocatable, device :: wa_dev(:,:)
   real, allocatable, device :: oa_dev(:,:)
   real, allocatable, device :: cwc_dev(:,:,:)
   real, allocatable, device :: fcld_dev(:,:)
   real, allocatable, device :: reff_dev(:,:,:)
   real, allocatable, device :: rsuvbm_dev(:)
   real, allocatable, device :: rsuvdf_dev(:)
   real, allocatable, device :: rsirbm_dev(:)
   real, allocatable, device :: rsirdf_dev(:)
   real, allocatable, device :: taua_dev(:,:,:)
   real, allocatable, device :: ssaa_dev(:,:,:)
   real, allocatable, device :: asya_dev(:,:,:)

   real, allocatable, device :: coa(:,:)
   real, allocatable, device :: cah(:,:)
   real, allocatable, device :: caib(:,:,:)
   real, allocatable, device :: caif(:,:)

   ! Device outputs 
   ! ---------------

   real, allocatable, device :: flx_dev(:,:)
   real, allocatable, device :: flc_dev(:,:)
   real, allocatable, device :: flxu_dev(:,:)
   real, allocatable, device :: flcu_dev(:,:)
   real, allocatable, device :: fdiruv_dev(:)
   real, allocatable, device :: fdifuv_dev(:)
   real, allocatable, device :: fdirpar_dev(:)
   real, allocatable, device :: fdifpar_dev(:)
   real, allocatable, device :: fdirir_dev(:)
   real, allocatable, device :: fdifir_dev(:)

   real, allocatable, device :: flx_sfc_band_dev(:,:)

   ! Constant Memory
   ! ---------------

   real, dimension(5), constant :: hk_uv
   real, dimension(5), constant :: wk_uv
   real, dimension(5), constant :: zk_uv
   real, dimension(5), constant :: ry_uv
   real, dimension(3), constant :: aig_uv
   real, dimension(3), constant :: awg_uv
   real, dimension(3), constant :: arg_uv
   real, constant :: aib_uv
   real, dimension(2), constant :: awb_uv
   real, dimension(2), constant :: arb_uv

   real, dimension(3,10), constant :: hk_ir
   real, dimension(10), constant :: xk_ir
   real, dimension(3), constant :: ry_ir
   real, constant :: aib_nir
   real, dimension(3,2), constant :: awb_nir
   real, dimension(3,2), constant :: arb_nir
   real, dimension(3,3), constant :: aia_nir
   real, dimension(3,3), constant :: awa_nir
   real, dimension(3,3), constant :: ara_nir
   real, dimension(3,3), constant :: aig_nir
   real, dimension(3,3), constant :: awg_nir
   real, dimension(3,3), constant :: arg_nir
#endif

   ! Parameters
   ! ----------

   integer, parameter :: nu = 43
   integer, parameter :: nw = 37
   integer, parameter :: nx = 62
   integer, parameter :: ny = 101
   integer, parameter :: nband_uv = 5
   integer, parameter :: nk_ir = 10
   integer, parameter :: nband_ir = 3

   integer, parameter :: nband = nband_uv + nband_ir

   real,    parameter :: dsm = 0.602

contains

#ifdef _CUDA
   attributes(global) subroutine sorad (m,np,nb,co2,ict,icb)
#else
   subroutine sorad (m,np,nb,cosz_dev,pl_dev,ta_dev,wa_dev,oa_dev,co2,&
         cwc_dev,fcld_dev,ict,icb,reff_dev,hk_uv,hk_ir,&
         taua_dev,ssaa_dev,asya_dev,&
         rsuvbm_dev,rsuvdf_dev,rsirbm_dev,rsirdf_dev,&
         flx_dev,flc_dev,fdiruv_dev,fdifuv_dev,&
         fdirpar_dev,fdifpar_dev,fdirir_dev,fdifir_dev,&
         flxu_dev,flcu_dev,&
         flx_sfc_band_dev&
)
#endif


!**********************************************************************
! Changes in November 2011:
!    Code rewritten to pass in the full taua, ssaa, and asya arrays
!    in echo of the old, original sorad.f. This is done to allow for
!    RRTMG aerosol work to be done similarly in the GC.

! Changes in July 2010:
!   
!   (1) Code relooped to have one, main loop over the soundings, many
!       arrays reduced in dimensionality, some to scalars.
!   (2) Instead of calling to external aerosol Mie code, efficiency and
!       asymmetry tables are now passed in.
! 
! Changes in August 2004:
!
!   (1) A layer was added above pl(1) to account for the absorption 
!       in the region 0-pl(1) hPa. pl(1) can be either 0 or >0 hPa.
!   (2) All parameters in the subroutine deledd were assigned to real*8
!   (3) The minimun values of water vapor and ozone, as well as tausto,
!       were changed to avoid computer precision problem. 
!
!***********************************************************************
!
!
! The maximum-random assumption is applied for treating cloud
!  overlapping. Clouds are grouped into high, middle, and low clouds
!  separated by the level indices ict and icb.  Note: ict must be 
! less than icb, and icb must be less than np+1.
!
! In a high spatial-resolution atmospheric model, fractional cloud cover
!  might be computed to be either 0 or 1.  In such a case, scaling of the
!  cloud optical thickness is not necessary, and the computation can be
!  made faster by compiling with -DOVERCAST.  Otherwise, do not compile
!  with the -DOVERCAST flag. (note: for the case that fractional cloud 
!  cover in a layer is either 0 or 1, the results of using either option
!  are identical).
!
!----- Input parameters (NOTE: _dev suffix is left off of input arrays
!                        for clarity's sake. _dev is needed to prevent
!                        GPU nameclash):                           
!                                                   units      size
!  number of soundings (m)                          n/d         1
!  number of atmospheric layers (np)                n/d         1
!  number of bands (nb)                             n/d         1   
!  cosine of solar zenith angle (cosz)              n/d         m
!  level pressure (pl)                              mb        m*(np+1)
!        pl(np+1) is the surface pressure
!  layer temperature (ta)                           k         m*np
!  layer specific humidity (wa)                     gm/gm     m*np
!  layer ozone concentration (oa)                   gm/gm     m*np
!  co2 mixing ratio by volume (co2)                 pppv        1
!  option for scaling cloud optical thickness       n/d         1
!        overcast="true" if scaling is NOT required
!        (the case that cloud cover is either 0 or 1)
!        overcast="false" if scaling is required
!  cloud water mixing ratio (cwc)                  gm/gm      m*np*3
!        index 1 for ice particle
!        index 2 for liquid drops
!        index 3 for rain drops
!        index 4 for snow
!  cloud amount (fcld)                             fraction   m*np
!  level index separating high and middle           n/d         1
!        clouds (ict)
!  level index separating middle and low            n/d         1
!        clouds (icb)
!  effective size of cloud particles (reff)        micron     m*np*3
!           index  1 for ice
!           index  2 for water drops
!           index  3 rain (not used in this code)
!           index  4 snow
!  aerosol optical thickness (taua)                 n/d      m*np*nb
!  aerosol single scattering albedo (ssaa)          n/d      m*np*nb
!  aerosol asymmetry factor (asya)                  n/d      m*np*nb
!  surface reflectivity     
!        in the UV+par region:
!           for beam insolation    (rsuvbm)        fraction     m 
!           for diffuse insolation (rsuvdf)        fraction     m 
!        in the near-ir region:
!           for beam insolation    (rsirbm)        fraction     m  
!           for diffuse insolation (rsirdf)        fraction     m
!
!  The 8 bands are:
!
!        in the uv region :
!           index  1 for the 0.225-0.285 micron band
!           index  2 for the 0.175-0.225;0.285-0.300 micron band
!           index  3 for the 0.300-0.325 micron band
!           index  4 for the 0.325-0.4 micron band
!        in the par region :
!           index  5 for the 0.4-0.690 micron band
!        in the infrared region :
!           index  6 for the 0.690-1.220 micron band
!           index  7 for the 1.220-2.270 micron band
!           index  8 for the 2.270-3.850 micron band
!
!----- Output parameters
!
!   all-sky flux (downward minus upward) (flx)     fraction   m*(np+1)
!   clear-sky flux (downward minus upward) (flc)   fraction   m*(np+1)
!   all-sky direct downward uv (0.175-0.4 micron)
!                flux at the surface (fdiruv)      fraction   m
!   all-sky diffuse downward uv flux at
!                the surface (fdifuv)              fraction   m
!   all-sky direct downward par (0.4-0.69 micron)
!                flux at the surface (fdirpar)     fraction   m
!   all-sky diffuse downward par flux at
!                the surface (fdifpar)             fraction   m
!   all-sky direct downward ir (0.69-10 micron)
!                flux at the surface (fdirir)      fraction   m
!   all-sky diffuse downward ir flux at
!                the surface (fdifir)              fraction   m
!   all-sky upward flux (flxu)                     fraction   m*(np+1)
!   clear-sky upward flux (flcu)                   fraction   m*(np+1)
!   all-sky flux (downward minus upward)           fraction   m*nb
!                per band at the surface
!                (flx_sfc_band)
!
!----- Notes:
!
!  (1) The unit of output fluxes (flx,flc,etc.) is fraction of the total
!      insolation at the top of the atmosphere.  Therefore, fluxes
!      are the output fluxes multiplied by the extra-terrestrial solar
!      flux and the cosine of the solar zenith angle.
!  (2) pl(*,1) is the pressure at the top of the model, and
!      pl(*,np+1) is the surface pressure.
!  (3) the pressure levels ict and icb correspond approximately
!      to 400 and 700 mb.
!
!-----Please notify Ming-Dah Chou for coding errors.
!        
!*************************************************************************
      implicit none


!-----input values

#ifdef _CUDA
      integer, value :: m,np,ict,icb,nb
      real, value :: co2
#else
!-----input parameters

      integer m,np,ict,icb,nb
      real cosz_dev(m),pl_dev(m,np+1),ta_dev(m,np),wa_dev(m,np),oa_dev(m,np),co2
      real cwc_dev(m,np,4),fcld_dev(m,np),reff_dev(m,np,4), hk_uv(5),hk_ir(3,10)
      real rsuvbm_dev(m),rsuvdf_dev(m),rsirbm_dev(m),rsirdf_dev(m)

      real taua_dev(m,np,nb)
      real ssaa_dev(m,np,nb)
      real asya_dev(m,np,nb)

!-----output parameters

      real flx_dev(m,np+1),flc_dev(m,np+1)
      real flxu_dev(m,np+1),flcu_dev(m,np+1)
      real fdiruv_dev (m),fdifuv_dev (m)
      real fdirpar_dev(m),fdifpar_dev(m)
      real fdirir_dev (m),fdifir_dev (m)
      real flx_sfc_band_dev(m,nband)
#endif

!-----temporary arrays

!GPU This is needed to allocate space for temporary local arrays on the GPU.
!GPU On the CPU this ifndef converts it back to np.

#ifndef GPU_MAXLEVS
#define GPU_MAXLEVS np
#endif

      integer :: i,k,l,ntop

      real :: dp(GPU_MAXLEVS),wh(GPU_MAXLEVS),oh(GPU_MAXLEVS)
      real :: scal(GPU_MAXLEVS)
      real :: swh(GPU_MAXLEVS+1),so2(GPU_MAXLEVS+1),df(0:GPU_MAXLEVS+1)
      real :: scal0, wvtoa, o3toa, pa
      real :: snt,cnt,xx4,xtoa
      real :: dp_pa(GPU_MAXLEVS)

!-----parameters for co2 transmission tables

      real :: w1,dw,u1,du

      integer :: ib
      real :: tauclb(GPU_MAXLEVS),tauclf(GPU_MAXLEVS),asycl(GPU_MAXLEVS)
      real :: taubeam(GPU_MAXLEVS,4),taudiff(GPU_MAXLEVS,4)
      real :: fcld_col(GPU_MAXLEVS)
      real :: cwc_col(GPU_MAXLEVS,4)
      real :: reff_col(GPU_MAXLEVS,4)
      real :: taurs,tauoz,tauwv
      real :: tausto,ssatau,asysto
      real :: tautob,ssatob,asytob
      real :: tautof,ssatof,asytof
      real :: rr(0:GPU_MAXLEVS+1,2),tt(0:GPU_MAXLEVS+1,2),td(0:GPU_MAXLEVS+1,2)
      real :: rs(0:GPU_MAXLEVS+1,2),ts(0:GPU_MAXLEVS+1,2)
      real :: fall(GPU_MAXLEVS+1),fclr(GPU_MAXLEVS+1),fsdir,fsdif
      real :: fupa(GPU_MAXLEVS+1),fupc(GPU_MAXLEVS+1)
      real :: cc1,cc2,cc3
      real :: rrt,ttt,tdt,rst,tst

      integer :: iv,ik
      real :: ssacl(GPU_MAXLEVS)

      integer :: im

      integer :: ic,iw 
      real :: ulog,wlog,dc,dd,x0,x1,x2,y0,y1,y2,du2,dw2

      integer :: ih
#ifdef OVERCAST
      real :: rra(0:GPU_MAXLEVS+1),rxa(0:GPU_MAXLEVS+1)
      real :: ttaold,tdaold,rsaold
      real :: ttanew,tdanew,rsanew
#else
      real :: rra(0:GPU_MAXLEVS+1,2,2),tta(0:GPU_MAXLEVS,2,2)
      real :: tda(0:GPU_MAXLEVS,2,2)
      real :: rsa(0:GPU_MAXLEVS,2,2),rxa(0:GPU_MAXLEVS+1,2,2)
#endif
      real :: flxdn
      real :: fdndir,fdndif,fupdif
      real :: denm,yy

      integer :: is
      real :: ch,cm,ct

      integer :: foundtop
      real :: dftop

      real :: dum

#ifdef _CUDA
!GPU--Initialize i for thread indexing
      i = (blockidx%x - 1) * blockdim%x + threadidx%x

      RUN_LOOP: if (i <= m) then
#else
      RUN_LOOP: do i=1,m
#endif

!GPU--GPU needs these variables initiated

         ntop = 0
         fdndir = 0.0
         fdndif = 0.0

!-----Beginning of sorad code

!-----wvtoa and o3toa are the water vapor and o3 amounts of the region 
!     above the pl(1) level.
!     snt is the secant of the solar zenith angle

         snt    = 1.0/cosz_dev(i)
         xtoa   = max(pl_dev(i,1),1.e-3)
         scal0  = xtoa*(0.5*xtoa/300.)**.8
         o3toa  = 1.02*oa_dev(i,1)*xtoa*466.7 + 1.0e-8
         wvtoa  = 1.02*wa_dev(i,1)*scal0 * (1.0+0.00135*(ta_dev(i,1)-240.)) + 1.0e-9
         swh(1)  = wvtoa

         do k=1,np

!-----compute layer thickness. indices for the surface level and
!     surface layer are np+1 and np, respectively.

            dp(k) = pl_dev(i,k+1)-pl_dev(i,k)
            dp_pa(k) = dp(k) * 100. ! dp in pascals
 
!-----compute scaled water vapor amount following Eqs. (3.3) and (3.5) 
!     unit is g/cm**2
!
            pa   = 0.5*(pl_dev(i,k)+pl_dev(i,k+1))
            scal(k) = dp(k)*(pa/300.)**.8
            wh(k)   = 1.02*wa_dev(i,k)*scal(k) * (1.+0.00135*(ta_dev(i,k)-240.)) + 1.e-9
            swh(k+1)= swh(k)+wh(k)

!-----compute ozone amount, unit is (cm-atm)stp
!     the number 466.7 is the unit conversion factor
!     from g/cm**2 to (cm-atm)stp

            oh(k)   = 1.02*oa_dev(i,k)*dp(k)*466.7 + 1.e-8

!-----Fill the reff, cwc, and fcld for the column

            fcld_col(k) = fcld_dev(i,k)
            do l = 1, 4
               reff_col(k,l) = reff_dev(i,k,l)
               cwc_col(k,l) = cwc_dev(i,k,l)
            end do

         end do

!-----Initialize temporary arrays to zero to avoid UMR

         rr = 0.0
         tt = 0.0
         td = 0.0
         rs = 0.0
         ts = 0.0

         rra = 0.0
         rxa = 0.0
#ifndef OVERCAST
         tta = 0.0
         tda = 0.0
         rsa = 0.0
#endif

!-----initialize fluxes for all-sky (flx), clear-sky (flc), and
!     flux reduction (df)
!
         do k=1,np+1
            flx_dev(i,k)=0.
            flc_dev(i,k)=0.
            flxu_dev(i,k)=0.
            flcu_dev(i,k)=0.
         end do

!-----Initialize new per-band surface fluxes

         do ib = 1, nband
            flx_sfc_band_dev(i,ib) = 0.
         end do

!-----Begin inline of SOLUV

!-----compute solar uv and par fluxes

!-----initialize fdiruv, fdifuv, surface reflectances and transmittances.
!     the reflectance and transmittance of the clear and cloudy portions
!     of a layer are denoted by 1 and 2, respectively.
!     cc is the maximum cloud cover in each of the high, middle, and low
!     cloud groups.
!     1/dsm=1/cos(53) = 1.66

         fdiruv_dev(i)=0.0
         fdifuv_dev(i)=0.0
         rr(np+1,1)=rsuvbm_dev(i)
         rr(np+1,2)=rsuvbm_dev(i)
         rs(np+1,1)=rsuvdf_dev(i)
         rs(np+1,2)=rsuvdf_dev(i)
         td(np+1,1)=0.0
         td(np+1,2)=0.0
         tt(np+1,1)=0.0
         tt(np+1,2)=0.0
         ts(np+1,1)=0.0
         ts(np+1,2)=0.0
         rr(0,1)=0.0
         rr(0,2)=0.0
         rs(0,1)=0.0
         rs(0,2)=0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
         tt(0,1)=1.0
         tt(0,2)=1.0
         ts(0,1)=1.0
         ts(0,2)=1.0
         cc1=0.0
         cc2=0.0
         cc3=0.0

!-----options for scaling cloud optical thickness

#ifdef OVERCAST

!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud asymmetry factor
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

         call getvistau(np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,0,0,&
                        taubeam,taudiff,asycl)

#else

!-----Compute scaled cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud asymmetry factor
!     Note: the cloud optical properties are assumed to be independent
!     of spectral bands in the UV and PAR regions.
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

         call getvistau(np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,ict,icb,&
                        taubeam,taudiff,asycl)

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group
!     The cc1,2,3 are still needed in the flux calculations below

!MAT---DO NOT FUSE THIS LOOP
!MAT---Loop must run to completion so that cc[1,2,3] are correct.
         do k=1,np
            if (k.lt.ict) then
               cc1=max(cc1,fcld_dev(i,k))
            elseif (k.lt.icb) then
               cc2=max(cc2,fcld_dev(i,k))
            else
               cc3=max(cc3,fcld_dev(i,k))
            end if
         end do
!MAT---DO NOT FUSE THIS LOOP

#endif

         do k=1,np
            tauclb(k)=taubeam(k,1)+taubeam(k,2)+taubeam(k,3)+taubeam(k,4)
            tauclf(k)=taudiff(k,1)+taudiff(k,2)+taudiff(k,3)+taudiff(k,4)
         end do

!-----integration over spectral bands

!-----Compute optical thickness, single-scattering albedo and asymmetry
!     factor for a mixture of "na" aerosol types. [Eqs. (4.16)-(4.18)]

         do ib=1,nband_uv

!-----compute direct beam transmittances of the layer above pl(1)

            td(0,1)=exp(-(wvtoa*wk_uv(ib)+o3toa*zk_uv(ib))/cosz_dev(i))
            td(0,2)=td(0,1)

            do k=1,np

!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor (Eqs. 6.2-6.4)

               taurs=ry_uv(ib)*dp(k)
               tauoz=zk_uv(ib)*oh(k)
               tauwv=wk_uv(ib)*wh(k)

               tausto=taurs+tauoz+tauwv+taua_dev(i,k,ib)+1.0e-7
               ssatau=ssaa_dev(i,k,ib)+taurs
               asysto=asya_dev(i,k,ib)

               tautob=tausto
               asytob=asysto/ssatau
               ssatob=ssatau/tautob+1.0e-8
               ssatob=min(ssatob,0.999999)

!-----for direct incident radiation

               call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

               call deledd(tautob,ssatob,asytob,dsm,rst,tst,dum)

               rr(k,1)=rrt
               tt(k,1)=ttt
               td(k,1)=tdt
               rs(k,1)=rst
               ts(k,1)=tst

!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer

!-----for direct incident radiation
!     The effective layer optical properties. Eqs. (6.2)-(6.4)

               tautob=tausto+tauclb(k)
               ssatob=(ssatau+tauclb(k))/tautob+1.0e-8
               ssatob=min(ssatob,0.999999)
               asytob=(asysto+asycl(k)*tauclb(k))/(ssatob*tautob)

!-----for diffuse incident radiation

               tautof=tausto+tauclf(k)
               ssatof=(ssatau+tauclf(k))/tautof+1.0e-8
               ssatof=min(ssatof,0.999999)
               asytof=(asysto+asycl(k)*tauclf(k))/(ssatof*tautof)

!-----for direct incident radiation
!     note that the cloud optical thickness is scaled differently 
!     for direct and diffuse insolation, Eqs. (7.3) and (7.4).

               call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation 
!     with an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

               call deledd(tautof,ssatof,asytof,dsm,rst,tst,dum)

               rr(k,2)=rrt
               tt(k,2)=ttt
               td(k,2)=tdt
               rs(k,2)=rst
               ts(k,2)=tst
            end do

!-----flux calculations
!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)

            do k=1,np+1
               fclr(k)=0.0
               fall(k)=0.0
               fupa(k)=0.0
               fupc(k)=0.0
            end do

            fsdir=0.0
            fsdif=0.0

#ifdef OVERCAST

!-----Inline CLDFLXY

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is either 0 or 1.

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated by 
!         beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)

!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

!-----ih=1 for clear sky; ih=2 for cloudy sky.

!-----First set is ih = 1
            rra(np+1)=rr(np+1,1)
            rxa(np+1)=rs(np+1,1)

            do k=np,0,-1
               denm=ts(k,1)/(1.-rs(k,1)*rxa(k+1))
               rra(k)=rr(k,1)+(td(k,1)*rra(k+1)+(tt(k,1)-td(k,1))*rxa(k+1))*denm
               rxa(k)=rs(k,1)+ts(k,1)*rxa(k+1)*denm
            end do

            tdaold = td(0,1)
            ttaold = tt(0,1)
            rsaold = rs(0,1)

            tdanew = 0.0
            ttanew = 0.0
            rsanew = 0.0

            do k=1,np+1
               if (k <= np) then
                  denm=ts(k,1)/(1.-rsaold*rs(k,1))
                  tdanew=tdaold*td(k,1)
                  ttanew=tdaold*tt(k,1)+(tdaold*rsaold*rr(k,1)+ttaold-tdaold)*denm
                  rsanew=rs(k,1)+ts(k,1)*rsaold*denm
               end if

               denm=1./(1.-rsaold*rxa(k))
               fdndir=tdaold
               xx4=tdaold*rra(k)
               yy=ttaold-tdaold
               fdndif=(xx4*rsaold+yy)*denm
               fupdif=(xx4+yy*rxa(k))*denm
               flxdn=fdndir+fdndif-fupdif
               fupc(k)=fupdif 
               fclr(k)=flxdn

               tdaold = tdanew
               ttaold = ttanew
               rsaold = rsanew

               tdanew = 0.0
               ttanew = 0.0
               rsanew = 0.0
            end do

!-----Second set is ih = 2

            rra(np+1)=rr(np+1,2)
            rxa(np+1)=rs(np+1,2)

            do k=np,0,-1
               denm=ts(k,2)/(1.-rs(k,2)*rxa(k+1))
               rra(k)=rr(k,2)+(td(k,2)*rra(k+1)+(tt(k,2)-td(k,2))*rxa(k+1))*denm
               rxa(k)=rs(k,2)+ts(k,2)*rxa(k+1)*denm
            end do

            tdaold = td(0,2)
            ttaold = tt(0,2)
            rsaold = rs(0,2)

            tdanew = 0.0
            ttanew = 0.0
            rsanew = 0.0

            do k=1,np+1
               if (k <= np) then
                  denm=ts(k,2)/(1.-rsaold*rs(k,2))
                  tdanew=tdaold*td(k,2)
                  ttanew=tdaold*tt(k,2)+(tdaold*rsaold*rr(k,2)+ttaold-tdaold)*denm
                  rsanew=rs(k,2)+ts(k,2)*rsaold*denm
               end if

               denm=1./(1.-rsaold*rxa(k))
               fdndir=tdaold
               xx4=tdaold*rra(k)
               yy=ttaold-tdaold
               fdndif=(xx4*rsaold+yy)*denm
               fupdif=(xx4+yy*rxa(k))*denm
               flxdn=fdndir+fdndif-fupdif

               fupa(k)=fupdif
               fall(k)=flxdn

               tdaold = tdanew
               ttaold = ttanew
               rsaold = rsanew

               tdanew = 0.0
               ttanew = 0.0
               rsanew = 0.0
            end do

            fsdir=fdndir
            fsdif=fdndif

!-----End CLDFLXY inline

#else

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is allowed to be between 0 and 1.
!     the all-sky flux, fall is the summation inside the brackets
!     of Eq. (7.11)

!-----Inline CLDFLX

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated 
!         by beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----To save memory space, tda, tta, and rsa are pre-computed 
!     for k<icb. The dimension of these parameters is (m,np,2,2). 
!     It would have been (m,np,2,2,2) if these parameters were 
!     computed for all k's.

!-----for high clouds
!     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition

            do ih=1,2
               tda(0,ih,1)=td(0,ih)
               tta(0,ih,1)=tt(0,ih)
               rsa(0,ih,1)=rs(0,ih)
               tda(0,ih,2)=td(0,ih)
               tta(0,ih,2)=tt(0,ih)
               rsa(0,ih,2)=rs(0,ih)

               do k=1,ict-1
                  denm=ts(k,ih)/(1.-rsa(k-1,ih,1)*rs(k,ih))
                  tda(k,ih,1)=tda(k-1,ih,1)*td(k,ih)
                  tta(k,ih,1)=tda(k-1,ih,1)*tt(k,ih)+(tda(k-1,ih,1)*rsa(k-1,ih,1)&
                        *rr(k,ih)+tta(k-1,ih,1)-tda(k-1,ih,1))*denm
                  rsa(k,ih,1)=rs(k,ih)+ts(k,ih)*rsa(k-1,ih,1)*denm
                  tda(k,ih,2)=tda(k,ih,1)
                  tta(k,ih,2)=tta(k,ih,1)
                  rsa(k,ih,2)=rsa(k,ih,1)
               end do ! k loop

!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition

               do k=ict,icb-1
                  do im=1,2
                     denm=ts(k,im)/(1.-rsa(k-1,ih,im)*rs(k,im))
                     tda(k,ih,im)=tda(k-1,ih,im)*td(k,im)
                     tta(k,ih,im)=tda(k-1,ih,im)*tt(k,im)+(tda(k-1,ih,im)&
                           *rsa(k-1,ih,im)*rr(k,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                     rsa(k,ih,im)=rs(k,im)+ts(k,im)*rsa(k-1,ih,im)*denm
                  end do ! im loop
               end do ! k loop
            end do ! ih loop

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----To save memory space, rra and rxa are pre-computed for k>=icb.
!     the dimension of these parameters is (m,np,2,2). It would have
!     been (m,np,2,2,2) if these parameters were computed for all k's.

!-----for the low clouds
!     is=1 for clear-sky condition, is=2 for cloudy-sky condition

            do is=1,2
               rra(np+1,1,is)=rr(np+1,is)
               rxa(np+1,1,is)=rs(np+1,is)
               rra(np+1,2,is)=rr(np+1,is)
               rxa(np+1,2,is)=rs(np+1,is)

               do k=np,icb,-1
                  denm=ts(k,is)/(1.-rs(k,is)*rxa(k+1,1,is))
                  rra(k,1,is)=rr(k,is)+(td(k,is)*rra(k+1,1,is)+(tt(k,is)-td(k,is))&
                        *rxa(k+1,1,is))*denm
                  rxa(k,1,is)=rs(k,is)+ts(k,is)*rxa(k+1,1,is)*denm
                  rra(k,2,is)=rra(k,1,is)
                  rxa(k,2,is)=rxa(k,1,is)
               end do ! k loop

!-----for middle clouds

               do k=icb-1,ict,-1
                  do im=1,2
                     denm=ts(k,im)/(1.-rs(k,im)*rxa(k+1,im,is))
                     rra(k,im,is)=rr(k,im)+(td(k,im)*rra(k+1,im,is)+(tt(k,im)-td(k,im))&
                           *rxa(k+1,im,is))*denm
                     rxa(k,im,is)=rs(k,im)+ts(k,im)*rxa(k+1,im,is)*denm
                  end do ! im loop
               end do ! k loop
            end do ! is loop

!-----integration over eight sky situations.
!     ih, im, is denote high, middle and low cloud groups.

            do ih=1,2

!-----clear portion 
               if(ih.eq.1) then
                  ch=1.0-cc1
!-----cloudy portion
               else
                  ch=cc1
               end if

               do im=1,2
!-----clear portion 
                  if(im.eq.1) then
                     cm=ch*(1.0-cc2)
!-----cloudy portion
                  else
                     cm=ch*cc2 
                  end if

                  do is=1,2
!-----clear portion 
                     if(is.eq.1) then
                        ct=cm*(1.0-cc3) 
!-----cloudy portion
                     else
                        ct=cm*cc3
                     end if

!-----add one layer at a time, going down.

                     do k=icb,np
                        denm=ts(k,is)/(1.-rsa(k-1,ih,im)*rs(k,is))
                        tda(k,ih,im)=tda(k-1,ih,im)*td(k,is)
                        tta(k,ih,im)=tda(k-1,ih,im)*tt(k,is)+(tda(k-1,ih,im)*rr(k,is)&
                              *rsa(k-1,ih,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                        rsa(k,ih,im)=rs(k,is)+ts(k,is)*rsa(k-1,ih,im)*denm
                     end do ! k loop

!-----add one layer at a time, going up.

                     do k=ict-1,0,-1
                        denm=ts(k,ih)/(1.-rs(k,ih)*rxa(k+1,im,is))
                        rra(k,im,is)=rr(k,ih)+(td(k,ih)*rra(k+1,im,is)+(tt(k,ih)-td(k,ih))&
                              *rxa(k+1,im,is))*denm
                        rxa(k,im,is)=rs(k,ih)+ts(k,ih)*rxa(k+1,im,is)*denm
                     end do ! k loop

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)

!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

                     do k=1,np+1
                        denm=1./(1.-rsa(k-1,ih,im)*rxa(k,im,is))
                        fdndir=tda(k-1,ih,im)
                        xx4=tda(k-1,ih,im)*rra(k,im,is)
                        yy=tta(k-1,ih,im)-tda(k-1,ih,im)
                        fdndif=(xx4*rsa(k-1,ih,im)+yy)*denm
                        fupdif=(xx4+yy*rxa(k,im,is))*denm
                        flxdn=fdndir+fdndif-fupdif

!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)

                        if(ih.eq.1 .and. im.eq.1 .and. is.eq.1) then
                           fupc(k)=fupdif
                           fclr(k)=flxdn
                        end if
                        fupa(k)=fupa(k)+fupdif*ct
                        fall(k)=fall(k)+flxdn*ct
                     end do ! k loop
                     fsdir=fsdir+fdndir*ct
                     fsdif=fsdif+fdndif*ct
                  end do ! is loop
               end do ! im loop
            end do ! ih loop

!-----End CLDFLX inline

#endif

!-----flux integration, Eq. (6.1)

            do k=1,np+1
               flx_dev(i,k)=flx_dev(i,k)+fall(k)*hk_uv(ib)
               flc_dev(i,k)=flc_dev(i,k)+fclr(k)*hk_uv(ib)
               flxu_dev(i,k)=flxu_dev(i,k)+fupa(k)*hk_uv(ib)
               flcu_dev(i,k)=flcu_dev(i,k)+fupc(k)*hk_uv(ib)
            end do

!-----get surface flux for each band
            flx_sfc_band_dev(i,ib)=flx_sfc_band_dev(i,ib)+fall(np+1)*hk_uv(ib)

!-----compute direct and diffuse downward surface fluxes in the UV
!     and par regions

            if(ib.lt.5) then
               fdiruv_dev(i)=fdiruv_dev(i)+fsdir*hk_uv(ib)
               fdifuv_dev(i)=fdifuv_dev(i)+fsdif*hk_uv(ib)
            else
               fdirpar_dev(i)=fsdir*hk_uv(ib)
               fdifpar_dev(i)=fsdif*hk_uv(ib)
            end if
         end do

!-----Inline SOLIR

!-----compute and update solar ir fluxes

         fdirir_dev(i)=0.0
         fdifir_dev(i)=0.0
         rr(np+1,1)=rsirbm_dev(i)
         rr(np+1,2)=rsirbm_dev(i)
         rs(np+1,1)=rsirdf_dev(i)
         rs(np+1,2)=rsirdf_dev(i)
         td(np+1,1)=0.0
         td(np+1,2)=0.0
         tt(np+1,1)=0.0
         tt(np+1,2)=0.0
         ts(np+1,1)=0.0
         ts(np+1,2)=0.0
         rr(0,1)=0.0
         rr(0,2)=0.0
         rs(0,1)=0.0
         rs(0,2)=0.0
!         td(0,1)=1.0
!         td(0,2)=1.0
         tt(0,1)=1.0
         tt(0,2)=1.0
         ts(0,1)=1.0
         ts(0,2)=1.0
         cc1=0.0
         cc2=0.0
         cc3=0.0

!-----integration over spectral bands

!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10)
!     The indices 1, 2, 3 are for ice, water, rain particles,
!     respectively.

         do ib=1,nband_ir
            iv=ib+5

!-----options for scaling cloud optical thickness

#ifdef OVERCAST

!-----Compute cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

            call getnirtau(ib,np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,0,0,&
                           taubeam,taudiff,asycl,ssacl)

#else

!-----Compute scaled cloud optical thickness. Eqs. (4.6) and (4.10) and
!     compute cloud single scattering albedo and asymmetry factor
!     for a mixture of ice and liquid particles.
!     Eqs.(4.6)-(4.8), (6.2)-(6.4)
!
!     The indices 1, 2, 3, 4 are for ice, water, rain, snow particles,
!     respectively.

            call getnirtau(ib,np,cosz_dev(i),dp_pa,fcld_col,reff_col,cwc_col,ict,icb,&
                           taubeam,taudiff,asycl,ssacl)

!-----clouds within each of the high, middle, and low clouds are 
!     assumed to be maximally overlapped, and the cloud cover (cc) 
!     for a group (high, middle, or low) is the maximum cloud cover 
!     of all the layers within a group

!MAT--DO NOT FUSE THIS LOOP
!MAT  Loop must run to completion so that cc[1,2,3] are correct.
            do k=1,np
               if (k.lt.ict) then
                  cc1=max(cc1,fcld_dev(i,k))
               elseif (k.lt.icb) then
                  cc2=max(cc2,fcld_dev(i,k))
               else
                  cc3=max(cc3,fcld_dev(i,k))
               end if
            end do
!MAT--DO NOT FUSE THIS LOOP

#endif

            do k=1,np
               tauclb(k)=taubeam(k,1)+taubeam(k,2)+taubeam(k,3)+taubeam(k,4)
               tauclf(k)=taudiff(k,1)+taudiff(k,2)+taudiff(k,3)+taudiff(k,4)
            end do


!-----integration over the k-distribution function

            do ik=1,nk_ir

!-----compute direct beam transmittances of the layer above pl(1)

               td(0,1)=exp(-wvtoa*xk_ir(ik)/cosz_dev(i))
               td(0,2)=td(0,1)

               do k=1,np
                  taurs=ry_ir(ib)*dp(k)
                  tauwv=xk_ir(ik)*wh(k)

!-----compute clear-sky optical thickness, single scattering albedo,
!     and asymmetry factor. Eqs.(6.2)-(6.4)

                  tausto=taurs+tauwv+taua_dev(i,k,iv)+1.0e-7
                  ssatau=ssaa_dev(i,k,iv)+taurs+1.0e-8
                  asysto=asya_dev(i,k,iv)
                  tautob=tausto
                  asytob=asysto/ssatau
                  ssatob=ssatau/tautob+1.0e-8
                  ssatob=min(ssatob,0.999999)

!-----Compute reflectance and transmittance of the clear portion 
!     of a layer

!-----for direct incident radiation

                  call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

                  call deledd(tautob,ssatob,asytob,dsm,rst,tst,dum)

                  rr(k,1)=rrt
                  tt(k,1)=ttt
                  td(k,1)=tdt
                  rs(k,1)=rst
                  ts(k,1)=tst

!-----compute reflectance and transmittance of the cloudy portion 
!     of a layer

!-----for direct incident radiation. Eqs.(6.2)-(6.4)

                  tautob=tausto+tauclb(k)
                  ssatob=(ssatau+ssacl(k)*tauclb(k))/tautob+1.0e-8
                  ssatob=min(ssatob,0.999999)
                  asytob=(asysto+asycl(k)*ssacl(k)*tauclb(k))/(ssatob*tautob)

!-----for diffuse incident radiation

                  tautof=tausto+tauclf(k)
                  ssatof=(ssatau+ssacl(k)*tauclf(k))/tautof+1.0e-8
                  ssatof=min(ssatof,0.999999)
                  asytof=(asysto+asycl(k)*ssacl(k)*tauclf(k))/(ssatof*tautof)

!-----for direct incident radiation

                  call deledd(tautob,ssatob,asytob,cosz_dev(i),rrt,ttt,tdt)

!-----diffuse incident radiation is approximated by beam radiation with
!     an incident angle of 53 degrees, Eqs.(6.5) and (6.6)

                  call deledd(tautof,ssatof,asytof,dsm,rst,tst,dum)

                  rr(k,2)=rrt
                  tt(k,2)=ttt
                  td(k,2)=tdt
                  rs(k,2)=rst
                  ts(k,2)=tst
               end do

!-----FLUX CALCULATIONS

!     initialize clear-sky flux (fclr), all-sky flux (fall), 
!     and surface downward fluxes (fsdir and fsdif)

               do k=1,np+1
                  fclr(k)=0.0
                  fall(k)=0.0
                  fupc(k)=0.0
                  fupa(k)=0.0
               end do

               fsdir=0.0
               fsdif=0.0

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is either 0 or 1.

#ifdef OVERCAST

!-----Inline CLDFLXY

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated by 
!         beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)

!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

!-----ih=1 for clear sky; ih=2 for cloudy sky.

!-----First set is ih = 1
               rra(np+1)=rr(np+1,1)
               rxa(np+1)=rs(np+1,1)

               do k=np,0,-1
                  denm=ts(k,1)/(1.-rs(k,1)*rxa(k+1))
                  rra(k)=rr(k,1)+(td(k,1)*rra(k+1)+(tt(k,1)-td(k,1))*rxa(k+1))*denm
                  rxa(k)=rs(k,1)+ts(k,1)*rxa(k+1)*denm
               end do

               tdaold = td(0,1)
               ttaold = tt(0,1)
               rsaold = rs(0,1)

               tdanew = 0.0
               ttanew = 0.0
               rsanew = 0.0

               do k=1,np+1
                  if (k <= np) then
                     denm=ts(k,1)/(1.-rsaold*rs(k,1))
                     tdanew=tdaold*td(k,1)
                     ttanew=tdaold*tt(k,1)+(tdaold*rsaold*rr(k,1)+ttaold-tdaold)*denm
                     rsanew=rs(k,1)+ts(k,1)*rsaold*denm
                  end if

                  denm=1./(1.-rsaold*rxa(k))
                  fdndir=tdaold
                  xx4=tdaold*rra(k)
                  yy=ttaold-tdaold
                  fdndif=(xx4*rsaold+yy)*denm
                  fupdif=(xx4+yy*rxa(k))*denm
                  flxdn=fdndir+fdndif-fupdif

                  fupc(k)=fupdif
                  fclr(k)=flxdn

                  tdaold = tdanew
                  ttaold = ttanew
                  rsaold = rsanew

                  tdanew = 0.0
                  ttanew = 0.0
                  rsanew = 0.0
               end do

!-----Second set is ih = 2

               rra(np+1)=rr(np+1,2)
               rxa(np+1)=rs(np+1,2)

               do k=np,0,-1
                  denm=ts(k,2)/(1.-rs(k,2)*rxa(k+1))
                  rra(k)=rr(k,2)+(td(k,2)*rra(k+1)+(tt(k,2)-td(k,2))*rxa(k+1))*denm
                  rxa(k)=rs(k,2)+ts(k,2)*rxa(k+1)*denm
               end do

               tdaold = td(0,2)
               ttaold = tt(0,2)
               rsaold = rs(0,2)

               tdanew = 0.0
               ttanew = 0.0
               rsanew = 0.0

               do k=1,np+1
                  if (k <= np) then
                     denm=ts(k,2)/(1.-rsaold*rs(k,2))
                     tdanew=tdaold*td(k,2)
                     ttanew=tdaold*tt(k,2)+(tdaold*rsaold*rr(k,2)+ttaold-tdaold)*denm
                     rsanew=rs(k,2)+ts(k,2)*rsaold*denm
                  end if

                  denm=1./(1.-rsaold*rxa(k))
                  fdndir=tdaold
                  xx4=tdaold*rra(k)
                  yy=ttaold-tdaold
                  fdndif=(xx4*rsaold+yy)*denm
                  fupdif=(xx4+yy*rxa(k))*denm
                  flxdn=fdndir+fdndif-fupdif

                  fupa(k)=fupdif
                  fall(k)=flxdn

                  tdaold = tdanew
                  ttaold = ttanew
                  rsaold = rsanew

                  tdanew = 0.0
                  ttanew = 0.0
                  rsanew = 0.0
               end do

               fsdir=fdndir
               fsdif=fdndif

!-----End CLDFLXY inline

#else

!-----for clear- and all-sky flux calculations when fractional 
!     cloud cover is allowed to be between 0 and 1.
!     the all-sky flux, fall is the summation inside the brackets
!     of Eq. (7.11)

!-----Inline CLDFLX

!-----compute transmittances and reflectances for a composite of
!     layers. layers are added one at a time, going down from the top.
!     tda is the composite direct transmittance illuminated 
!         by beam radiation
!     tta is the composite total transmittance illuminated by
!         beam radiation
!     rsa is the composite reflectance illuminated from below
!         by diffuse radiation
!     tta and rsa are computed from Eqs. (6.10) and (6.12)

!-----To save memory space, tda, tta, and rsa are pre-computed 
!     for k<icb. The dimension of these parameters is (m,np,2,2). 
!     It would have been (m,np,2,2,2) if these parameters were 
!     computed for all k's.

!-----for high clouds
!     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition

               do ih=1,2
                  tda(0,ih,1)=td(0,ih)
                  tta(0,ih,1)=tt(0,ih)
                  rsa(0,ih,1)=rs(0,ih)
                  tda(0,ih,2)=td(0,ih)
                  tta(0,ih,2)=tt(0,ih)
                  rsa(0,ih,2)=rs(0,ih)

                  do k=1,ict-1
                     denm=ts(k,ih)/(1.-rsa(k-1,ih,1)*rs(k,ih))
                     tda(k,ih,1)=tda(k-1,ih,1)*td(k,ih)
                     tta(k,ih,1)=tda(k-1,ih,1)*tt(k,ih)+(tda(k-1,ih,1)*rsa(k-1,ih,1)&
                           *rr(k,ih)+tta(k-1,ih,1)-tda(k-1,ih,1))*denm
                     rsa(k,ih,1)=rs(k,ih)+ts(k,ih)*rsa(k-1,ih,1)*denm
                     tda(k,ih,2)=tda(k,ih,1)
                     tta(k,ih,2)=tta(k,ih,1)
                     rsa(k,ih,2)=rsa(k,ih,1)
                  end do ! k loop

!-----for middle clouds
!     im=1 for clear-sky condition, im=2 for cloudy-sky condition

                  do k=ict,icb-1
                     do im=1,2
                        denm=ts(k,im)/(1.-rsa(k-1,ih,im)*rs(k,im))
                        tda(k,ih,im)=tda(k-1,ih,im)*td(k,im)
                        tta(k,ih,im)=tda(k-1,ih,im)*tt(k,im)+(tda(k-1,ih,im)&
                              *rsa(k-1,ih,im)*rr(k,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                        rsa(k,ih,im)=rs(k,im)+ts(k,im)*rsa(k-1,ih,im)*denm
                     end do ! im loop
                  end do ! k loop
               end do ! ih loop

!-----layers are added one at a time, going up from the surface.
!     rra is the composite reflectance illuminated by beam radiation
!     rxa is the composite reflectance illuminated from above
!         by diffuse radiation
!     rra and rxa are computed from Eqs. (6.9) and (6.11)

!-----To save memory space, rra and rxa are pre-computed for k>=icb.
!     the dimension of these parameters is (m,np,2,2). It would have
!     been (m,np,2,2,2) if these parameters were computed for all k's.

!-----for the low clouds
!     is=1 for clear-sky condition, is=2 for cloudy-sky condition

               do is=1,2
                  rra(np+1,1,is)=rr(np+1,is)
                  rxa(np+1,1,is)=rs(np+1,is)
                  rra(np+1,2,is)=rr(np+1,is)
                  rxa(np+1,2,is)=rs(np+1,is)

                  do k=np,icb,-1
                     denm=ts(k,is)/(1.-rs(k,is)*rxa(k+1,1,is))
                     rra(k,1,is)=rr(k,is)+(td(k,is)*rra(k+1,1,is)+(tt(k,is)-td(k,is))&
                           *rxa(k+1,1,is))*denm
                     rxa(k,1,is)=rs(k,is)+ts(k,is)*rxa(k+1,1,is)*denm
                     rra(k,2,is)=rra(k,1,is)
                     rxa(k,2,is)=rxa(k,1,is)
                  end do ! k loop

!-----for middle clouds

                  do k=icb-1,ict,-1
                     do im=1,2
                        denm=ts(k,im)/(1.-rs(k,im)*rxa(k+1,im,is))
                        rra(k,im,is)=rr(k,im)+(td(k,im)*rra(k+1,im,is)+(tt(k,im)-td(k,im))&
                              *rxa(k+1,im,is))*denm
                        rxa(k,im,is)=rs(k,im)+ts(k,im)*rxa(k+1,im,is)*denm
                     end do ! im loop
                  end do ! k loop
               end do ! is loop

!-----integration over eight sky situations.
!     ih, im, is denote high, middle and low cloud groups.

               do ih=1,2
!-----clear portion 
                  if(ih.eq.1) then
                     ch=1.0-cc1
!-----cloudy portion
                  else
                     ch=cc1
                  end if

                  do im=1,2
!-----clear portion 
                     if(im.eq.1) then
                        cm=ch*(1.0-cc2)
!-----cloudy portion
                     else
                        cm=ch*cc2 
                     end if

                     do is=1,2
!-----clear portion 
                        if(is.eq.1) then
                           ct=cm*(1.0-cc3) 
!-----cloudy portion
                        else
                           ct=cm*cc3
                        end if

!-----add one layer at a time, going down.

                        do k=icb,np
                           denm=ts(k,is)/(1.-rsa(k-1,ih,im)*rs(k,is))
                           tda(k,ih,im)=tda(k-1,ih,im)*td(k,is)
                           tta(k,ih,im)=tda(k-1,ih,im)*tt(k,is)+(tda(k-1,ih,im)*rr(k,is)&
                                 *rsa(k-1,ih,im)+tta(k-1,ih,im)-tda(k-1,ih,im))*denm
                           rsa(k,ih,im)=rs(k,is)+ts(k,is)*rsa(k-1,ih,im)*denm
                        end do ! k loop

!-----add one layer at a time, going up.

                        do k=ict-1,0,-1
                           denm=ts(k,ih)/(1.-rs(k,ih)*rxa(k+1,im,is))
                           rra(k,im,is)=rr(k,ih)+(td(k,ih)*rra(k+1,im,is)+(tt(k,ih)-td(k,ih))&
                                 *rxa(k+1,im,is))*denm
                           rxa(k,im,is)=rs(k,ih)+ts(k,ih)*rxa(k+1,im,is)*denm
                        end do ! k loop

!-----compute fluxes following Eq. (6.15) for fupdif and
!     Eq. (6.16) for (fdndir+fdndif)
!     fdndir is the direct  downward flux
!     fdndif is the diffuse downward flux
!     fupdif is the diffuse upward flux

                        do k=1,np+1
                           denm=1./(1.-rsa(k-1,ih,im)*rxa(k,im,is))
                           fdndir=tda(k-1,ih,im)
                           xx4=tda(k-1,ih,im)*rra(k,im,is)
                           yy=tta(k-1,ih,im)-tda(k-1,ih,im)
                           fdndif=(xx4*rsa(k-1,ih,im)+yy)*denm
                           fupdif=(xx4+yy*rxa(k,im,is))*denm
                           flxdn=fdndir+fdndif-fupdif

!-----summation of fluxes over all sky situations;
!     the term in the brackets of Eq. (7.11)

                           if(ih.eq.1 .and. im.eq.1 .and. is.eq.1) then
                              fupc(k)=fupdif
                              fclr(k)=flxdn
                           end if
                           fupa(k)=fupa(k)+fupdif*ct
                           fall(k)=fall(k)+flxdn*ct
                        end do ! k loop
                        fsdir=fsdir+fdndir*ct
                        fsdif=fsdif+fdndif*ct
                     end do ! is loop
                  end do ! im loop
               end do ! ih loop

!-----End CLDFLX inline
#endif

!-----flux integration following Eq. (6.1)

               do k=1,np+1
                  flx_dev(i,k)=flx_dev(i,k)+fall(k)*hk_ir(ib,ik)
                  flc_dev(i,k)=flc_dev(i,k)+fclr(k)*hk_ir(ib,ik)
                  flxu_dev(i,k)=flxu_dev(i,k)+fupa(k)*hk_ir(ib,ik)
                  flcu_dev(i,k)=flcu_dev(i,k)+fupc(k)*hk_ir(ib,ik)
               end do

!-----compute downward surface fluxes in the ir region

               fdirir_dev(i)=fdirir_dev(i)+fsdir*hk_ir(ib,ik)
               fdifir_dev(i)=fdifir_dev(i)+fsdif*hk_ir(ib,ik)

!-----tabulate surface flux at ir bands
               flx_sfc_band_dev(i,iv)=flx_sfc_band_dev(i,iv)+fall(np+1)*hk_ir(ib,ik)

            end do ! ik loop
         end do

!-----compute pressure-scaled o2 amount following Eq. (3.5) with f=1.
!     unit is (cm-atm)stp. 165.22 = (1000/980)*23.14%*(22400/32)
!     compute flux reduction due to oxygen following Eq. (3.18). 0.0633 is the
!     fraction of insolation contained in the oxygen bands

         df(0) = 0.0
         cnt = 165.22*snt
         so2(1) = scal0*cnt
! LLT increased parameter 145 to 155 to enhance effect
         df(1) = 0.0633*(1.-exp(-0.000155*sqrt(so2(1))))

         do k=1,np
            so2(k+1) = so2(k) + scal(k)*cnt
! LLT increased parameter 145 to 155 to enhance effect
            df(k+1) = 0.0633*(1.0 - exp(-0.000155*sqrt(so2(k+1)))) 
         end do

!-----for solar heating due to co2 scaling follows Eq(3.5) with f=1.
!     unit is (cm-atm)stp. 789 = (1000/980)*(44/28.97)*(22400/44)

         so2(1) = (789.*co2)*scal0

         do k=1,np
            so2(k+1) = so2(k) + (789.*co2)*scal(k)
         end do

!-----The updated flux reduction for co2 absorption in Band 7 where absorption due to
!     water vapor and co2 are both moderate. df is given by the second term on the
!     right-hand-side of Eq. (3.24) divided by So. so2 and swh are the co2 and
!     water vapor amounts integrated from the top of the atmosphere

         u1 = -3.0
         du = 0.15
         w1 = -4.0
         dw = 0.15

!-----Inline RFLX
         du2=du*du
         dw2=dw*dw
         x0=u1+real(nu)*du
         y0=w1+real(nw)*dw

         x1=u1-0.5*du
         y1=w1-0.5*dw

         do k= 1, np+1
            ulog=min(log10(so2(k)*snt),x0)
            wlog=min(log10(swh(k)*snt),y0)
            ic=int((ulog-x1)/du+1.)
            iw=int((wlog-y1)/dw+1.)
            if(ic.lt.2)  ic=2
            if(iw.lt.2)  iw=2
            if(ic.gt.nu) ic=nu
            if(iw.gt.nw) iw=nw
            dc=ulog-real(ic-2)*du-u1
            dd=wlog-real(iw-2)*dw-w1   
            x2=cah(ic-1,iw-1)+(cah(ic-1,iw)-cah(ic-1,iw-1))/dw*dd
            y2=x2+(cah(ic,iw-1)-cah(ic-1,iw-1))/du*dc
            y2=max(y2,0.0)
            df(k)=df(k)+1.5*y2 ! LLT increase CO2 effect to help reduce cold tropopause bias
         end do

!-----df is the updated flux reduction for co2 absorption
!     in Band 8 where the co2 absorption has a large impact
!     on the heating of middle atmosphere. From the table
!     given by Eq. (3.19)

         u1 = 0.000250
         du = 0.000050
         w1 = -2.0
         dw = 0.05

!-----Inline RFLX
         du2=du*du
         dw2=dw*dw
         x0=u1+real(nx)*du
         y0=w1+real(ny)*dw

         x1=u1-0.5*du
         y1=w1-0.5*dw

         do k= 1, np+1
            ulog=min(co2*snt,x0)
            wlog=min(log10(pl_dev(i,k)),y0)
            ic=int((ulog-x1)/du+1.)
            iw=int((wlog-y1)/dw+1.)
            if(ic.lt.2)  ic=2
            if(iw.lt.2)  iw=2
            if(ic.gt.nx) ic=nx
            if(iw.gt.ny) iw=ny
            dc=ulog-real(ic-2)*du-u1
            dd=wlog-real(iw-2)*dw-w1   
            x2=coa(ic-1,iw-1)+(coa(ic-1,iw)-coa(ic-1,iw-1))/dw*dd
            y2=x2+(coa(ic,iw-1)-coa(ic-1,iw-1))/du*dc
            y2=max(y2,0.0)
            df(k)=df(k)+1.5*y2 ! LLT increase CO2 effect to help reduce cold tropopause bias
         end do

!-----adjust the o2-co2 reduction below cloud top following Eq. (6.18)

!GPU--A new routine for finding ntop and using it
!GPU  is provided. This helps shield the GPUs from the
!GPU  exit statement which causes divergent warps.

         foundtop = 0

         do k=1,np
            if (fcld_dev(i,k) > 0.02.and.foundtop.eq.0) then
               foundtop = 1
               ntop = k
            end if
         end do

         if (foundtop.eq.0) ntop=np+1

         dftop = df(ntop)

         do k=1,np+1
            if (k .gt. ntop) then
               xx4   = (flx_dev(i,k)/flx_dev(i,ntop))
               df(k) = dftop + xx4 * (df(k)-dftop)
            end if
         end do

!-----update the net fluxes

         do k=1,np+1
            df(k) = min(df(k),flx_dev(i,k)-1.0e-8)
!           df(k) = 0.0
            flx_dev(i,k) = flx_dev(i,k) - df(k)
            flc_dev(i,k) = flc_dev(i,k) - df(k)
         end do

!-----update the downward surface fluxes 

!        xx4 = fdirir (i) + fdifir (i) +&
!              fdiruv (i) + fdifuv (i) +&
!              fdirpar(i) + fdifpar(i)

         xx4 = flx_dev(i,np+1) + df(np+1)

         if ( abs(xx4) > epsilon(1.0) ) then
            xx4 = max(min(1.0 - df(np+1)/xx4,1.),0.)
         else
            xx4 = 0.0
         end if

         fdirir_dev(i)  = xx4*fdirir_dev(i) 
         fdifir_dev(i)  = xx4*fdifir_dev(i) 
         fdiruv_dev(i)  = xx4*fdiruv_dev(i) 
         fdifuv_dev(i)  = xx4*fdifuv_dev(i) 
         fdirpar_dev(i) = xx4*fdirpar_dev(i)
         fdifpar_dev(i) = xx4*fdifpar_dev(i)

         do ib = 1, nband
            flx_sfc_band_dev(i,ib) = xx4*flx_sfc_band_dev(i,ib)
         end do

#ifndef _CUDA
      end do RUN_LOOP
#else
      end if RUN_LOOP
#endif

   end subroutine sorad

!*********************************************************************

#ifdef _CUDA
   attributes(device) subroutine deledd(tau1,ssc1,g01,cza1,rr1,tt1,td1)
#else
   subroutine deledd(tau1,ssc1,g01,cza1,rr1,tt1,td1)
#endif

!*********************************************************************
!
!-----uses the delta-eddington approximation to compute the
!     bulk scattering properties of a single layer
!     coded following King and Harshvardhan (JAS, 1986)
!
!  inputs:
!     tau1:  optical thickness
!     ssc1:  single scattering albedo
!     g01:   asymmetry factor
!     cza1:  cosine o the zenith angle
!
!  outputs:
!
!     rr1:  reflection of the direct beam
!     tt1:  total (direct+diffuse) transmission of the direct beam
!     td1:  direct transmission of the direct beam
!
!*********************************************************************

      implicit none

!-----input parameters

      real(kind=MAPL_R4), intent(in) :: tau1,ssc1,g01,cza1

!-----output parameters

      real(kind=MAPL_R4), intent(out) :: rr1, tt1, td1

!-----temporary parameters

      real(kind=MAPL_R8), parameter :: zero = 0.0_MAPL_R8
      real(kind=MAPL_R8), parameter :: one = 1.0_MAPL_R8
      real(kind=MAPL_R8), parameter :: two = 2.0_MAPL_R8
      real(kind=MAPL_R8), parameter :: three = 3.0_MAPL_R8
      real(kind=MAPL_R8), parameter :: four = 4.0_MAPL_R8
      real(kind=MAPL_R8), parameter :: fourth = 0.25_MAPL_R8
      real(kind=MAPL_R8), parameter :: seven = 7.0_MAPL_R8
      real(kind=MAPL_R8), parameter :: thresh = 1.e-8_MAPL_R8

      real(kind=MAPL_R8) ::  tau,ssc,g0,rr,tt,td 
      real(kind=MAPL_R8) ::  zth,ff,xx,taup,sscp,gp,gm1,gm2,gm3,akk,alf1,alf2
      real(kind=MAPL_R8) ::  all,bll,st7,st8,cll,dll,fll,ell,st1,st2,st3,st4

      zth = real(cza1,kind=MAPL_R8)
      g0  = real(g01 ,kind=MAPL_R8)
      tau = real(tau1,kind=MAPL_R8)
      ssc = real(ssc1,kind=MAPL_R8)

      ff  = g0*g0
      xx  = one-ff*ssc
      taup= tau*xx
      sscp= ssc*(one-ff)/xx
      gp  = g0/(one+g0)

      xx  = three*gp 
      gm1 = (seven-sscp*(four+xx))*fourth
      gm2 =-(one  -sscp*(four-xx))*fourth

      akk = sqrt((gm1+gm2)*(gm1-gm2))

      xx  = akk*zth
      st7 = one-xx
      st8 = one+xx
      st3 = st7*st8

      if (abs(st3) .lt. thresh) then
         zth = zth+0.0010
         if(zth > 1.0) zth = zth-0.0020 
         xx  = akk*zth
         st7 = one-xx
         st8 = one+xx
         st3 = st7*st8
      end if

      td=exp(-taup/zth)

      gm3 = (two-zth*three*gp)*fourth
      xx  = gm1-gm2
      alf1= gm1-gm3*xx
      alf2= gm2+gm3*xx

      xx  = akk*two
      all = (gm3-alf2*zth    )*xx*td
      bll = (one-gm3+alf1*zth)*xx

      xx  = akk*gm3
      cll = (alf2+xx)*st7
      dll = (alf2-xx)*st8

      xx  = akk*(one-gm3)
      fll = (alf1+xx)*st8
      ell = (alf1-xx)*st7

      st2 = exp(-akk*taup)
      st4 = st2*st2

      st1 = sscp/((akk+gm1+(akk-gm1)*st4)*st3)

      rr = ( cll-dll*st4     - all*st2)*st1
      tt =-((fll-ell*st4)*td - bll*st2)*st1

      rr = max(rr,zero)
      tt = max(tt,zero)

      tt = tt+td

      td1 = real(td,kind=MAPL_R4)
      rr1 = real(rr,kind=MAPL_R4)
      tt1 = real(tt,kind=MAPL_R4)

   end subroutine deledd

! GPU Due to limitations of CUDA Fortran, we cannot call device subroutines
!     external to the file containing the global subroutine. So we replicate
!     the interface for the tau routine here.

#ifdef _CUDA
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: getvistau
!
! !INTERFACE:
   attributes(device) subroutine getvistau(nlevs,cosz,dp,fcld,reff,hydromets,ict,icb,&
                                           taubeam,taudiff,asycl)

      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: nlevs          !  Number of levels
      real,    intent(IN ) :: cosz           !  Cosine of solar zenith angle
      real,    intent(IN ) :: dp(:)          !  Delta pressure (Pa)
      real,    intent(IN ) :: fcld(:)        !  Cloud fraction (used sometimes)
      real,    intent(IN ) :: reff(:,:)      !  Effective radius (microns)
      real,    intent(IN ) :: hydromets(:,:) !  Hydrometeors (kg/kg)
      integer, intent(IN ) :: ict, icb       !  Flags for various uses 
!                 ict  = 0   Indicates that in-cloud values have been given
!                            and are expected
!                 ict != 0   Indicates that overlap computation is needed, and:
!                               ict is the level of the mid-high boundary
!                               icb is the level of the low-mid  boundary
!                
! !OUTPUT PARAMETERS:
      real,    intent(OUT) :: taubeam(:,:)   !  Optical Depth for Beam Radiation
      real,    intent(OUT) :: taudiff(:,:)   !  Optical Depth for Diffuse Radiation
      real,    intent(OUT) ::   asycl(:  )   !  Cloud Asymmetry Factor
! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for visible wavelengths
!  In general will compute in-cloud - will do grid mean when called
!  for diagnostic use only. ict flag will indicate which to do.
!  Slots for reff, hydrometeors, taubeam, taudiff, and asycl are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.10.27   Molod moved to Radiation_Shared and revised arg list, units
!    2011.11.16   MAT: Generalized to a call that is per-column
!
!EOP
!------------------------------------------------------------------------------
!BOC

      integer            :: k,in,im,it,ia,kk
      real               :: fm,ft,fa,xai,tauc,asyclt
      real               :: cc(3)
      real               :: taucld1,taucld2,taucld3,taucld4
      real               :: g1,g2,g3,g4

      real               :: reff_snow

#include "getvistau.code"

      return

!EOC
   end subroutine getvistau
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: getnirtau
!
! !INTERFACE:
   attributes(device) subroutine getnirtau(ib,nlevs,cosz,dp,fcld,reff,hydromets,ict,icb,&
                                           taubeam,taudiff,asycl,ssacl)

      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: ib             !  Band number
      integer, intent(IN ) :: nlevs          !  Number of levels
      real,    intent(IN ) :: cosz           !  Cosine of solar zenith angle
      real,    intent(IN ) :: dp(:)          !  Delta pressure in Pa
      real,    intent(IN ) :: fcld(:)        !  Cloud fraction (used sometimes)
      real,    intent(IN ) :: reff(:,:)      !  Effective radius (microns)
      real,    intent(IN ) :: hydromets(:,:) !  Hydrometeors (kg/kg)
      integer, intent(IN ) :: ict, icb       !  Flags for various uses 
!                 ict  = 0   Indicates that in-cloud values have been given
!                            and are expected
!                 ict != 0   Indicates that overlap computation is needed, and:
!                               ict is the level of the mid-high boundary
!                               icb is the level of the low-mid  boundary
!                
! !OUTPUT PARAMETERS:
      real,    intent(OUT) :: taubeam(:,:)   !  Optical depth for beam radiation
      real,    intent(OUT) :: taudiff(:,:)   !  Optical depth for diffuse radiation
      real,    intent(OUT) ::   ssacl(:  )   !  Cloud single scattering albedo
      real,    intent(OUT) ::   asycl(:  )   !  Cloud asymmetry factor
! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for near-infrared wavelengths
!  In general will compute in-cloud - will do grid mean when called
!  for diagnostic use only. ict flag will indicate which to do.
!  Slots for reff, hydrometeors and tauall are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.10.27   Molod moved to Radiation_Shared and revised arg list, units
!    2011.11.16   MAT: Generalized to a call that is per-column
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer            :: k,in,im,it,ia,kk
      real               :: fm,ft,fa,xai,tauc,asyclt,ssaclt
      real               :: cc(3)
      real               :: taucld1,taucld2,taucld3,taucld4
      real               :: g1,g2,g3,g4
      real               :: w1,w2,w3,w4

      real               :: reff_snow

#include "getnirtau.code"

      return

!EOC
   end subroutine getnirtau

#endif

end module
