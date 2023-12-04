module LockEntrain
! <OVERVIEW>
!
!      K-PROFILE BOUNDARY LAYER SCHEME WITH CLOUD TOP ENTRAINMENT
!      Adapted from code provided by Stephen A. Klein at GFDL
!
!      This routine calculates diffusivity coefficients for vertical
!      diffusion using a K-profile approach.  This scheme is modelled
!      after:
!
!      Lock, A.P., A.R. Brown, M.R. Bush, G.M. Martin, and R.N.B. Smith, 
!          2000: A new boundary layer mixing scheme. Part I: Scheme 
!          description and single-column modeling tests. Mon. Wea. Rev.,
!          128, 3187-3199.
!
!   
! </OVERVIEW>
! <DESCRIPTION>
!
!      The key part is the parameterization of entrainment at the top
!      convective layers. For an entrainment interface from surface
!      driven mixing, the entrainment rate, we, is parameterized as:
!
!                      
!      we, surf =  A / B
!
!      where A = ( beta_surf * (V_surf**3 + V_shear**3) / zsml )
!        and B = ( delta_b   + ((V_surf**3 + V_shear**3)**(2/3))/zsml )
!
!
!      In this formula,
!
!           zsml     =  depth of surface mixed layer
!
!           V_surf   =  surface driven scaling velocity
!                    =  (u_star*b_star*zsml)**(1/3)
!
!           V_shear  =  surface driven shear velocity,
!                    =  (Ashear**(1/3))*u_star
!
!           delta_b  =  buoyancy jump at the entrainment interface(m/s2)
!                    =  grav * delta_slv / slv
!
!      If an entrainment interface is associated only with cloud top
!      radiative cooling, the entrainment rate is parameterized as:
!
!
!                     
!      we, rad  =  ( A / B)
!
!            where A = beta_rad  *  V_rad**3 /  zradml
!              and B = delta_b   +  V_rad**2 /  zradml
!
!      where
!
!           zradml   =  depth of radiatively driven layer
!
!           V_rad    =  radiatively driven scaling velocity
!                    =  (grav*delta-F*zradml/(rho*cp_air*T)) **(1/3)
!                 
!
!      Note that the source of entrainment from cloud top buoyancy
!      reversal has been omitted in this implementation.
!  
!      If the entrainment interface for surface driven mixing coincides
!      with that for cloud top radiatively driven convection then the
!      following full entrainment rate:
!
!                          
!      we, full =   A / B
!
!            where A =   V_full**3 / zsml
!              and B =  delta_b+((V_surf**3+V_shear**3+V_rad**3)**(2/3))/zsml
!              and V_full**3 = beta_surf*(V_surf**3+V_shear**3) + beta_rad*V_rad**3
!   
! </DESCRIPTION>
!

!-----------------------------------------------------------------------
!
! outside modules used 
!

#ifndef _CUDA
   use GEOS_UtilsMod, only: DQSAT=>GEOS_DQSAT
#else
   use cudafor
   ! NOTE: GPUs use the QSAT and DQSAT at the end of this module
#endif 

   use MAPL_ConstantsMod, only: MAPL_GRAV,  MAPL_KARMAN, MAPL_CP,     &
                                MAPL_RGAS,  MAPL_RVAP,   MAPL_ALHL,   &
                                MAPL_ALHS,  MAPL_TICE,   MAPL_VIREPS, &
                                MAPL_P00,   MAPL_KAPPA,  MAPL_H2OMW,  &
                                MAPL_AIRMW, MAPL_R4,     MAPL_R8
   use MAPL,              only: MAPL_UNDEF

   implicit none

#ifndef _CUDA
   private

!-----------------------------------------------------------------------
!
!  public interfaces

   PUBLIC ENTRAIN

#endif

#ifdef _CUDA

   ! Inputs
   ! ------

   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: TDTLW_IN_dev
   REAL, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: U_STAR_dev
   REAL, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: B_STAR_dev
   REAL, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: FRLAND_dev
   REAL, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: EVAP_dev
   REAL, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: SH_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: T_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: QV_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: QL_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: QI_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: U_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: V_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: ZFULL_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: PFULL_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: ZHALF_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: PHALF_dev

   ! Inoutputs
   ! ---------

   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: DIFF_M_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: DIFF_T_dev

   ! Outputs
   ! -------

   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: K_M_ENTR_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: K_T_ENTR_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: K_SFC_dev
   REAL, ALLOCATABLE, DIMENSION(:,:,:), DEVICE :: K_RAD_dev
   REAL, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: ZCLOUD_dev
   REAL, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: ZRADML_dev
   REAL, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: ZRADBASE_dev
   REAL, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: ZSML_dev

   ! Diagnostics
   ! -----------

   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: ZCLDTOP_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: WENTR_SFC_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: WENTR_RAD_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: DEL_BUOY_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: VSFC_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: VRAD_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: KENTRAD_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: VBRV_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: WENTR_BRV_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: DSIEMS_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: CHIS_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: DELSINV_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: SLMIXTURE_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: CLDRADF_DIAG_dev
   REAL, TARGET, ALLOCATABLE, DIMENSION(:,:  ), DEVICE :: RADRCODE_DIAG_dev

   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: ZCLDTOP_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: WENTR_SFC_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: WENTR_RAD_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: DEL_BUOY_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: VSFC_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: VRAD_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: KENTRAD_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: VBRV_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: WENTR_BRV_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: DSIEMS_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: CHIS_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: DELSINV_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: SLMIXTURE_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: CLDRADF_DIAG_dev_ptr
   REAL, POINTER,             DIMENSION(:,:  ), DEVICE :: RADRCODE_DIAG_dev_ptr

   ! Parameters for Internal DQSAT
   ! -----------------------------

   REAL, PARAMETER :: ESFAC            = MAPL_H2OMW/MAPL_AIRMW
   REAL, PARAMETER :: MAX_MIXING_RATIO = 1. 
   REAL, PARAMETER :: ZEROC            = MAPL_TICE

   REAL, PARAMETER :: TMINTBL   =  150.0
   REAL, PARAMETER :: TMAXTBL   =  333.0
   REAL, PARAMETER :: DEGSUBS   =  100  
   REAL, PARAMETER :: ERFAC     = (DEGSUBS/ESFAC)
   REAL, PARAMETER :: DELTA_T   =  1.0 / DEGSUBS
   REAL, PARAMETER :: TABLESIZE =  NINT(TMAXTBL-TMINTBL)*DEGSUBS + 1
   REAL, PARAMETER :: TMIX      = -20. 

   REAL, PARAMETER :: TMINSTR = -95. 
   REAL, PARAMETER :: TSTARR1 = -75. 
   REAL, PARAMETER :: TSTARR2 = -65. 
   REAL, PARAMETER :: TSTARR3 = -50. 
   REAL, PARAMETER :: TSTARR4 = -40. 
   REAL, PARAMETER :: TMAXSTR = +60. 

   REAL(KIND=MAPL_R8), PARAMETER :: B6 = 6.136820929E-11*100.0
   REAL(KIND=MAPL_R8), PARAMETER :: B5 = 2.034080948E-8 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: B4 = 3.031240396E-6 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: B3 = 2.650648471E-4 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: B2 = 1.428945805E-2 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: B1 = 4.436518521E-1 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: B0 = 6.107799961E+0 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: BI6= 1.838826904E-10*100.0
   REAL(KIND=MAPL_R8), PARAMETER :: BI5= 4.838803174E-8 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: BI4= 5.824720280E-6 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: BI3= 4.176223716E-4 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: BI2= 1.886013408E-2 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: BI1= 5.034698970E-1 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: BI0= 6.109177956E+0 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S16= 0.516000335E-11*100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S15= 0.276961083E-8 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S14= 0.623439266E-6 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S13= 0.754129933E-4 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S12= 0.517609116E-2 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S11= 0.191372282E+0 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S10= 0.298152339E+1 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S26= 0.314296723E-10*100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S25= 0.132243858E-7 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S24= 0.236279781E-5 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S23= 0.230325039E-3 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S22= 0.129690326E-1 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S21= 0.401390832E+0 *100.0
   REAL(KIND=MAPL_R8), PARAMETER :: S20= 0.535098336E+1 *100.0

   REAL(KIND=MAPL_R8), PARAMETER :: DI(0:3) = (/ 57518.5606E08, 2.01889049, 3.56654, 20.947031 /)
   REAL(KIND=MAPL_R8), PARAMETER :: CI(0:3) = (/ 9.550426, -5723.265, 3.53068, -.00728332 /)
   REAL(KIND=MAPL_R8), PARAMETER :: DL(1:6) = (/ -7.902980, 5.02808, -1.3816, 11.344, 8.1328, -3.49149 /)
   REAL(KIND=MAPL_R8), PARAMETER :: LOGPS   = 3.005714898  ! LOG10(1013.246)
   REAL(KIND=MAPL_R8), PARAMETER :: TS      = 373.16
   REAL(KIND=MAPL_R8), PARAMETER :: CL(0:9) = (/54.842763, -6763.22, -4.21000, .000367, &
                                       .0415, 218.8,  53.878000, -1331.22, -9.44523, .014025  /)

   REAL, PARAMETER :: TMINLQU = MAPL_TICE - 40.0
   REAL, PARAMETER :: TMINICE = MAPL_TICE + -95.

#endif

!-----------------------------------------------------------------------
!
!  set default values to some instance independent parameters       
!

   real, parameter :: akmax      =  1.e4 ! maximum value for a diffusion coefficient 
                                         ! (m2/s)

   real, parameter :: zcldtopmax =  3.e3 ! maximum altitude for cloud top of 
                                         ! radiatively driven convection (m)    

   real, parameter :: ramp       =  20.


!-----------------------------------------------------------------------
!
! declare version number 
!

   character(len=128) :: Version = '$Id$'
   character(len=128) :: Tagname = '$Name$'

!-----------------------------------------------------------------------
!
! Subroutines include:
!
!     entrain          main driver program of the module
!
!     mpbl_depth       routine to calculate the depth of surface driven
!                      mixed layer
!
!     radml_depth      subroutine to calculate the depth of the cloud
!                      topped radiatively driven mixed layer
!     
!     diffusivity_pbl2 subroutine to calculate diffusivity coefficients
!                      for surface driven mixed layer

contains

!======================================================================= 

#ifdef _CUDA
   attributes(global)    &
#endif 
   subroutine entrain(   &
! Integer Inputs
         icol,           &
         jcol,           &
         nlev,           &
! Inputs
         tdtlw_in,       &
         u_star,         &
         b_star,         &
         frland,         &
         evap,           &
         sh,             &
         t,              &
         qv,             &
         qlls,           &
         qils,           &
         u,              &
         v,              &
         zfull,          &
         pfull,          &
         zhalf,          &
         phalf,          &
! Inoutputs
         diff_m,         &
         diff_t,         &
! Outputs
         k_m_entr,       &
         k_t_entr,       &
         k_sfc,          &
         k_rad,          &
         zcloud,         &
         zradml,         &
         zradbase,       &
         zsml,           &
! Diagnostics
         zcldtop_diag,   &
         wentr_sfc_diag, &
         wentr_rad_diag, &
         del_buoy_diag,  &
         vsfc_diag,      &
         vrad_diag,      &
         kentrad_diag,   &
         vbrv_diag,      &
         wentr_brv_diag, &
         dsiems_diag,    &
         chis_diag,      &
         delsinv_diag,   &
         slmixture_diag, &
         cldradf_diag,   &
         radrcode_diag,  &
! Constants
         prandtlsfc,     &
         prandtlrad,     &
         beta_surf,      &
         beta_rad,       &
         tpfac_sfc,      &
         entrate_sfc,    &
         pceff_sfc,      &
         vscale_sfc,     &
         pertopt_sfc,    &
         khradfac,       &
         khsfcfac_lnd,   &
         khsfcfac_ocn    )

!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
!
!      ncol,nlev  sizes of i,j,k dimensions
!
!      convect   is surface based moist convection occurring in this
!                grid box?
!      u_star    friction velocity (m/s)
!      b_star    buoyancy scale (m/s**2)
!
!      three dimensional fields on model full levels, reals dimensioned
!      (:,:,pressure), third index running from top of atmosphere to 
!      bottom
!          
!      t         temperature (K)
!      qv        water vapor specific humidity (kg vapor/kg air)
!      ql        liquid water specific humidity (kg cond/kg air)
!      qi        ice water specific humidity (kg cond/kg air)
!      zfull     height of full levels (m)
!      pfull     pressure (Pa)
!      u         zonal wind (m/s)
!      v         meridional wind (m/s)
!
!      the following two fields are on the model half levels, with
!      size(zhalf,3) = size(t,3) +1, zhalf(:,:,size(zhalf,3)) 
!      must be height of surface (if you are not using eta-model)
!
!      zhalf     height at half levels (m)
!      phalf     pressure at half levels (Pa)
!
!      ------------
!      input/output
!      ------------
!
!      the following variables are defined at half levels and are
!      dimensions 1:nlev
!
!      diff_t   input and output heat diffusivity (m2/sec)
!      diff_m   input and output momentum diffusivity (m2/sec)
!
!      The diffusivity coefficient output from the routine includes
!      the modifications to use the internally calculated diffusivity
!      coefficients.
!
!      ------
!      output
!      ------
!
!      The following variables are defined at half levels and are
!      dimensions 1:nlev.
!
!      k_t_entr  heat diffusivity coefficient (m**2/s)
!      k_m_entr  momentum diffusivity coefficient (m**2/s)
!      zsml      height of surface driven mixed layer (m)
!
!      --------
!      internal
!      --------
!
!
!      General variables
!      -----------------
!
!      slv         virtual static energy (J/kg)      
!      density     air density (kg/m3)
!      hleff       effective latent heat of vaporization/sublimation 
!                  (J/kg)
!
!
!
!      Variables related to surface driven convective layers
!      -----------------------------------------------------
!
!      vsurf       surface driven buoyancy velocity scale (m/s)
!      vshear      surface driven shear velocity scale (m/s)
!      wentr_pbl   surface driven entrainment rate (m/s)
!      convpbl     1 is surface driven convective layer present
!                  0 otherwise
!      pblfq       1 if the half level is part of a surface driven
!                  layer, 0 otherwise
!      k_m_troen   momentum diffusion coefficient (m2/s) !Why are there two
!      k_t_troen   heat diffusion coefficient (m2/s)     !of these. They are equal?
!
!
!      Variables related to cloud top driven radiatively driven layers
!      ---------------------------------------------------------------
!
!      zradbase    height of base of radiatively driven mixed layer (m)
!      zradtop     height of top of radiatively driven mixed layer (m)
!      zradml      depth of radiatively driven mixed layer (m)
!      vrad        radiatively driven velocity scale (m/s)
!      radf        longwave jump at cloud top (W/m2) -- the radiative 
!                  forcing for cloud top driven mixing.
!      wentr_rad   cloud top driven entrainment (m/s)
!      svpcp       cloud top value of liquid water virtual static energy 
!                  divided by cp (K)
!      radpbl      1 if cloud top radiatively driven layer is present
!                  0 otherwise
!      radfq       1 if the half level is part of a radiatively driven
!                  layer, 0 otherwise
!      k_rad       radiatively driven diffusion coefficient (m2/s)
!
!-----------------------------------------------------------------------

      implicit none

#ifdef _CUDA
      integer, value,  intent(in) :: icol,jcol,nlev

      real,    value,  intent(in) :: prandtlsfc,prandtlrad,beta_surf,beta_rad
      real,    value,  intent(in) :: khradfac,tpfac_sfc,entrate_sfc,vscale_sfc,pertopt_sfc
      real,    value,  intent(in) :: pceff_sfc,khsfcfac_lnd,khsfcfac_ocn

      real,    device, intent(in),    dimension(icol,jcol,nlev)      :: tdtlw_in       
      real,    device, intent(in),    dimension(icol,jcol)           :: u_star,b_star,frland,evap,sh
      real,    device, intent(in),    dimension(icol,jcol,nlev)      :: t,qv,qlls,qils
      real,    device, intent(in),    dimension(icol,jcol,nlev)      :: u,v,zfull,pfull
      real,    device, intent(in),    dimension(icol,jcol,1:nlev+1)  :: zhalf, phalf ! 0:72 in GC, 1:73 here.
      real,    device, intent(inout), dimension(icol,jcol,1:nlev+1)  :: diff_m,diff_t
      real,    device, intent(out),   dimension(icol,jcol,1:nlev+1)  :: k_m_entr,k_t_entr
      real,    device, intent(out),   dimension(icol,jcol,1:nlev+1)  :: k_rad,k_sfc
      real,    device, intent(out),   dimension(icol,jcol)           :: zsml,zradml,zcloud,zradbase

      real,    device, pointer, dimension(:,:) :: wentr_rad_diag, wentr_sfc_diag ,del_buoy_diag
      real,    device, pointer, dimension(:,:) :: vrad_diag, kentrad_diag,vbrv_diag,wentr_brv_diag
      real,    device, pointer, dimension(:,:) :: dsiems_diag, chis_diag,zcldtop_diag
      real,    device, pointer, dimension(:,:) :: delsinv_diag, slmixture_diag,cldradf_diag
      real,    device, pointer, dimension(:,:) :: vsfc_diag,radrcode_diag
#else
      integer, intent(in)                                    :: icol,jcol,nlev

      real,    intent(in),    dimension(icol,jcol,nlev)      :: tdtlw_in       
      real,    intent(in),    dimension(icol,jcol)           :: u_star,b_star,frland,evap,sh
      real,    intent(in),    dimension(icol,jcol,nlev)      :: t,qv,qlls,qils
      real,    intent(in),    dimension(icol,jcol,nlev)      :: u,v,zfull,pfull
      real,    intent(in),    dimension(icol,jcol,1:nlev+1)  :: zhalf, phalf ! 0:72 in GC, 1:73 here.
      real,    intent(inout), dimension(icol,jcol,1:nlev+1)  :: diff_m,diff_t
      real,    intent(out),   dimension(icol,jcol,1:nlev+1)  :: k_m_entr,k_t_entr
      real,    intent(out),   dimension(icol,jcol,1:nlev+1)  :: k_rad,k_sfc
      real,    intent(out),   dimension(icol,jcol)           :: zsml,zradml,zcloud,zradbase

      real,    intent(in) :: prandtlsfc,prandtlrad,beta_surf,beta_rad
      real,    intent(in) :: khradfac,tpfac_sfc,entrate_sfc, vscale_sfc, pertopt_sfc
      real,    intent(in) :: pceff_sfc,khsfcfac_lnd,khsfcfac_ocn

      real, pointer, dimension(:,:) :: wentr_rad_diag, wentr_sfc_diag ,del_buoy_diag
      real, pointer, dimension(:,:) :: vrad_diag, kentrad_diag,vbrv_diag,wentr_brv_diag
      real, pointer, dimension(:,:) :: dsiems_diag, chis_diag,zcldtop_diag
      real, pointer, dimension(:,:) :: delsinv_diag, slmixture_diag,cldradf_diag
      real, pointer, dimension(:,:) :: vsfc_diag,radrcode_diag
#endif

! GPU The GPUs need to know how big local arrays are during compile-time
!     as the GPUs cannot allocate memory themselves. This command resets
!     this a priori size to nlev for the CPU.
#ifndef GPU_MAXLEVS
#define GPU_MAXLEVS nlev
#endif

      integer                      :: i,j,k,ibot
      integer                      :: ipbl
      integer                      :: kmax,kcldtop,kcldbot,kcldtop2
      logical                      :: do_jump_exit
      real                         :: maxradf,stab
      real                         :: wentrmax
      real                         :: qlcrit,ql00,qlm1,Abuoy,Ashear
      real                         :: wentr_tmp,hlf,vbulkshr,vbulk_scale
      real                         :: k_entr_tmp,tmpjump,critjump,radperturb,buoypert
      real                         :: tmp1, tmp2, slmixture
      real                         :: vsurf3, vshear3,vrad3, vbr3,dsiems
      real                         :: ztmp
      real                         :: zradtop
      real                         :: vrad,radf,svpcp,chis,vbrv
      real                         :: vsurf
      real                         :: wentr_rad,wentr_pbl,wentr_brv
      real                         :: convpbl, radpbl
      real, dimension(GPU_MAXLEVS) :: slv, density,qc, slvcp
      !real, dimension(GPU_MAXLEVS) :: rh ! Not used in current code
      real, dimension(GPU_MAXLEVS) :: hleff
      real, dimension(GPU_MAXLEVS) :: radfq,pblfq
      real, dimension(GPU_MAXLEVS) :: k_m_troen,k_t_troen
      real, dimension(GPU_MAXLEVS) :: k_rad_col
      real, dimension(GPU_MAXLEVS) :: dqs
      !real, dimension(GPU_MAXLEVS) :: qs ! Not used in current code
      real                         :: qs_dummy

!-----------------------------------------------------------------------
!
!     initialize variables

!-----------------------------------------------------------------------
!
!     Sizes
!

      qlcrit     = 1.0e-6
      Abuoy      = 0.23
      Ashear     = 25.0
      wentrmax   = 0.05

      ibot = nlev

#ifdef _CUDA
      i = (blockidx%x - 1) * blockdim%x + threadidx%x
      j = (blockidx%y - 1) * blockdim%y + threadidx%y

      I_LOOP: if ( i <= icol ) then
      J_LOOP: if ( j <= jcol ) then
#else
      I_LOOP: do i = 1, icol
      J_LOOP: do j = 1, jcol
#endif
         zradml(i,j)   = 0.0
         zcloud(i,j)   = 0.0
         zradbase(i,j) = 0.0
         zsml(i,j)     = 0.0 ! note that this must be zero as this is 
                             ! indicates stable surface layer and this
                             ! value is output for use in gravity
                             ! wave drag scheme
         zradtop    = 0.0
         convpbl    = 0.0
         wentr_pbl  = 0.0
         vsurf      = 0.0
         radpbl     = 0.0
         svpcp      = 0.0
         vrad       = 0.0
         radf       = 0.0
         wentr_rad  = 0.0
         K_LOOP_0: do k = 1, nlev
            pblfq(k)     = 0.0
            k_t_troen(k) = 0.0
            k_m_troen(k) = 0.0
            radfq(k)     = 0.0
            k_rad_col(k) = 0.0
            k_t_entr(i,j,k)   = 0.0
            k_m_entr(i,j,k)   = 0.0
            k_sfc(i,j,k) = 0.0
            k_rad(i,j,k) = 0.0
         end do K_LOOP_0

         ! For GPU we must zero out even the LM+1'th position

         k_t_entr(i,j,nlev+1)   = 0.0
         k_m_entr(i,j,nlev+1)   = 0.0
         k_sfc(i,j,nlev+1) = 0.0
         k_rad(i,j,nlev+1) = 0.0

!--------------------------------------------------------------------------
! Initialize optional outputs

         if(associated(wentr_sfc_diag)) wentr_sfc_diag(i,j) = MAPL_UNDEF
         if(associated(wentr_rad_diag)) wentr_rad_diag(i,j) = MAPL_UNDEF
         if(associated(del_buoy_diag))  del_buoy_diag(i,j)  = MAPL_UNDEF
         if(associated(vrad_diag))      vrad_diag(i,j)      = MAPL_UNDEF
         if(associated(vsfc_diag))      vsfc_diag(i,j)      = MAPL_UNDEF
         if(associated(kentrad_diag))   kentrad_diag(i,j)   = MAPL_UNDEF  
         if(associated(chis_diag))      chis_diag(i,j)      = MAPL_UNDEF
         if(associated(vbrv_diag))      vbrv_diag(i,j)      = MAPL_UNDEF
         if(associated(dsiems_diag))    dsiems_diag(i,j)    = MAPL_UNDEF
         if(associated(wentr_brv_diag)) wentr_brv_diag(i,j) = MAPL_UNDEF
         if(associated(zcldtop_diag))   zcldtop_diag(i,j)   = MAPL_UNDEF
         if(associated(slmixture_diag)) slmixture_diag(i,j) = MAPL_UNDEF
         if(associated(delsinv_diag))   delsinv_diag(i,j)   = MAPL_UNDEF
         if(associated(cldradf_diag))   cldradf_diag(i,j)   = MAPL_UNDEF
         if(associated(radrcode_diag))  radrcode_diag(i,j)  = MAPL_UNDEF

!-----------------------------------------------------------------------
!
!     set up specific humidities and static energies  
!     compute airdensity
!

         do k = 1, nlev
            if ( t(i,j,k) <= MAPL_TICE-ramp ) then
               HLEFF(k) = MAPL_ALHS
            else if ( (t(i,j,k) > MAPL_TICE-ramp) .and. (t(i,j,k) < MAPL_TICE) ) then
               HLEFF(k) =  ( (t(i,j,k)-MAPL_TICE+ramp)*MAPL_ALHL + &
                             (MAPL_TICE -t(i,j,k)    )*MAPL_ALHS   ) / ramp
            else
               HLEFF(k) = MAPL_ALHL
            end if

!--------------------------------------------------------------------------
!     Compute: 
!      qs        saturation specific humidity (kg/kg)
!      dqs       derivative of qs w/ respect to temp 

            ! Note. The array QS is never needed in the current code. So we
            ! use a dummy. QS can be restored if needed.
            !dqs(k) = dqsat(t(i,j,k), pfull(i,j,k), qsat=qs(k), pascals=.true. )
            dqs(k) = dqsat(t(i,j,k), pfull(i,j,k), qsat=qs_dummy, pascals=.true. )

!--------------------------------------------------------------------------
!     Compute total cloud condensate - qc. These are grid box mean values.

            qc(k) = (qlls(i,j,k) + qils(i,j,k))

!--------------------------------------------------------------------------
!     Compute relative humidity

            ! rh is not used in the current code. Commented out below.
            !rh(k) = qv(i,j,k) / qs(k)

!--------------------------------------------------------------------------
!     Compute liquid static energy.

            slv(k) = MAPL_CP*t(i,j,k)*(1+MAPL_VIREPS*qv(i,j,k)-qc(k)) + &
                         MAPL_GRAV*zfull(i,j,k) - hleff(k)*qc(k)

!!       slv(k)     = MAPL_CP*t(i,j,k) + MAPL_GRAV*zfull(i,j,k) - hleff(k)*qc(k)
!!       slv(k)     = slv(k)*(1+MAPL_VIREPS*(qv(i,j,k)+qc(k)))

            density(k) = pfull(i,j,k)/(MAPL_RGAS*(t(i,j,k) *(1.+MAPL_VIREPS*qv(i,j,k)-qc(k))))              

         end do

!--------------------------
! 
!     big loop over points
!

!---------------
! reset indices

         ipbl    = -1
         kcldtop = -1

!-----------------------------------------------------------
!
! SURFACE DRIVEN CONVECTIVE LAYERS
!
! Note this part is done only if b_star > 0., that is,
! upward surface buoyancy flux


         IF_BSTAR_GT_0: if (b_star(i,j) .gt. 0.) then

!------------------------------------------------------
! Find depth of surface driven mixing by raising a 
! parcel from the surface with some excess buoyancy
! to its level of neutral buoyancy.  Note the use
! slv as the density variable permits one to goes
! through phase changes to find parcel top

            call mpbl_depth(i,j,icol,jcol,nlev,&
                  tpfac_sfc,        &
                  entrate_sfc,      &
                  pceff_sfc,        &
                  vscale_sfc,       & 
                  pertopt_sfc,      &
                  t,                &
                  qv,               &
                  u,                &
                  v,                &
                  zfull,            &
                  pfull,            &
                  b_star,           &
                  u_star,           &
                  evap,             &
                  sh,               &
                  ipbl,zsml         )

!------------------------------------------------------
! Define velocity scales vsurf and vshear
!           
! vsurf   =  (u_star*b_star*zsml)**(1/3)
! vshear  =  (Ashear**(1/3))*u_star

            vsurf3   = u_star(i,j)*b_star(i,j)*zsml(i,j)
            vshear3  = Ashear*u_star(i,j)*u_star(i,j)*u_star(i,j)

            vsurf = vsurf3**(1./3.)
            if(associated(vsfc_diag)) vsfc_diag(i,j) = vsurf

!------------------------------------------------------
! Define velocity scale vbulkshr
!           

            vbulkshr   =  sqrt ( ( u(i,j,ipbl) - u(i,j,ibot) )**2        &
                               + ( v(i,j,ipbl) - v(i,j,ibot) )**2  )     &
                               /  zsml(i,j)
            vbulk_scale = 3.0 / 1000.   ! Non-local (PBL-deep) shear scale


!------------------------------------------------------
! Following Lock et al. 2000, limit height of surface
! well mixed layer if interior stable interface is
! found.  An interior stable interface is diagnosed if
! the slope between 2 full levels is greater than critjump

            critjump = 2.0
            if (ipbl .lt. ibot) then 
               do k = ibot, ipbl+1, -1
                  tmpjump =(slv(k-1)-slv(k))/MAPL_CP 
                  if (tmpjump .gt. critjump) then
                     ipbl = k
                     zsml(i,j) = zhalf(i,j,ipbl)
                     exit
                  end if
               enddo
            end if

!-------------------------------------
! compute entrainment rate
!

            tmp1 = MAPL_GRAV*max(0.1,(slv(ipbl-1)-slv(ipbl))/ &
                                 MAPL_CP)/(slv(ipbl)/MAPL_CP)
            tmp2 = ((vsurf3+vshear3)**(2./3.)) / zsml(i,j)

            wentr_tmp= min( wentrmax,  max(0., (beta_surf *   &
                           (vsurf3 + vshear3)/zsml(i,j))/ &
                           (tmp1+tmp2) ) )

!----------------------------------------
! fudgey adjustment of entrainment to reduce it
! for shallow boundary layers, and increase for 
! deep ones
            if ( zsml(i,j) .lt. 1600. ) then 
               wentr_tmp = wentr_tmp * ( zsml(i,j) / 800. )
            else
               wentr_tmp = 2.*wentr_tmp
            endif
!-----------------------------------------

!!AMM106 !----------------------------------------
!!AMM106 ! More fudgey adjustment of entrainment.
!!AMM106 ! Zeroes entr if bulk shear in PBL > vbulk_scale
!!AMM106 if ( vbulkshr .gt. vbulk_scale ) wentr_tmp = 0.0
!!AMM106 if (    ( vbulkshr .gt. 0.5*vbulk_scale  )   &
!!AMM106   .and. ( vbulkshr .le.     vbulk_scale  ) ) then 
!!AMM106      wentr_tmp = wentr_tmp * ( vbulk_scale -  vbulkshr ) *2 &
!!AMM106                                  / vbulk_scale
!!AMM106 endif

            k_entr_tmp = wentr_tmp*(zfull(i,j,ipbl-1)-zfull(i,j,ipbl))  
            k_entr_tmp = min ( k_entr_tmp, akmax )

            do k = ipbl, ibot
               pblfq(k) = 1.
            end do
            convpbl                 = 1.
            wentr_pbl               = wentr_tmp
            k_t_troen(ipbl)         = k_entr_tmp
            k_m_troen(ipbl)         = k_entr_tmp
            k_t_entr(i,j,ipbl)  = k_t_entr(i,j,ipbl) + k_entr_tmp
            k_m_entr(i,j,ipbl)  = k_m_entr(i,j,ipbl) + k_entr_tmp

            if(associated(wentr_sfc_diag)) wentr_sfc_diag(i,j) = wentr_tmp

!------------------------------------------------------
! compute diffusion coefficients in the interior of
! the PBL

            if (ipbl .lt. ibot) then

               call diffusivity_pbl2(i,j,icol,jcol,nlev, &
                     zsml,                  &
                     khsfcfac_lnd, khsfcfac_ocn, k_entr_tmp,  & 
                     vsurf,frland,          & 
                     zhalf,                 &
                     k_m_troen,             &
                     k_t_troen              )

               do k = ipbl+1, ibot
                  k_t_entr(i,j,k) =       & 
                        k_t_entr(i,j,k) + &
                        k_t_troen(k)
                        
                  k_m_entr(i,j,k) =       & 
                        k_m_entr(i,j,k) + &
                        k_m_troen(k)*prandtlsfc
               end do

            end if

         end if IF_BSTAR_GT_0


!-----------------------------------------------------------
!
! NEGATIVELY BUOYANT PLUMES DRIVEN BY 
! LW RADIATIVE COOLING AND/OR BUOYANCY REVERSAL
!
! This part is done only if a level kcldtop can be 
! found with: 
!
!    qc(kcldtop)>=qlcrit.and.qc(kcldtop-1)<qlcrit
!
! below zcldtopmax

         kmax = ibot+1
         do k = 1, ibot
            if( zhalf(i,j,k) < zcldtopmax) then
               kmax = k
               exit
            end if
         end do

!-----------------------------------------------------------
! Find cloud top and bottom using GRID BOX MEAN or IN-CLOUD 
! value of qc.  Decision occurs where qc is calculated


         kcldtop  = ibot+1
         do k = ibot,kmax,-1
            qlm1 = qc(k-1)  ! qc one level UP
            ql00 = qc(k)
            stab = slv(k-1) - slv(k) 
            if ( ( ql00  .ge. qlcrit ) .and. ( qlm1 .lt. qlcrit) .and. (stab .gt. 0.) ) then
               kcldtop  = k   
               exit
            end if
         end do

         if (kcldtop .ge. ibot+1) then 
            if(associated(radrcode_diag)) radrcode_diag(i,j)=1.
            go to 55
         endif

         kcldtop2=min( kcldtop+1,nlev)
         ! Look one level further down in case first guess is a thin diffusive veil
         if( (qc(kcldtop) .lt. 10*qlcrit ) .and. (qc(kcldtop2) .ge. 10*qc(kcldtop) ) ) then
            kcldtop=kcldtop2
         endif


         kcldbot  = ibot+1
         do k = ibot,kcldtop,-1
            qlm1 = qc(k-1)  ! qc one level UP
            ql00 = qc(k)
            if ( ( ql00  .lt. qlcrit ) .and. ( qlm1 .ge. qlcrit) ) then
               kcldbot  = k   
               exit
            end if
         end do


         if (kcldtop .eq. kcldbot) then 
            if(associated(radrcode_diag)) radrcode_diag(i,j)=2.
            go to 55 
         endif


         ! With diffusion of ql, qi "cloud top" found via these quantities may be above radiation max
         kcldtop2=min( kcldtop+2,nlev)
         maxradf = maxval( -1.*tdtlw_in(i,j,kcldtop:kcldtop2) )

         maxradf = maxradf*MAPL_CP*( (phalf(i,j,kcldtop+1)-phalf(i,j,kcldtop)) / MAPL_GRAV )

         maxradf = max( maxradf , 0. ) ! do not consider cloud tops that are heating

!!AMM108 ! Break out of "rad" plume if layer above cloud is 
!!AMM108 ! wetter than RH = somenumber
!!AMM108 if ( rh(kcldtop-1) .gt. 0.5 ) then 
!!AMM108    if(associated(radrcode_diag)) if radrcode_diag(i,j)=5.
!!AMM108    go to 55 
!!AMM108 endif  


!-----------------------------------------------------------
! Calculate optimal mixing fraction - chis - for buoyancy 
! reversal.  Use effective heat of evap/subl *tion.  Ignore 
! diffs across cldtop
         hlf = hleff(kcldtop)

         tmp1 = ( slv(kcldtop-1)  -  hlf*qc(kcldtop-1) ) - &
                ( slv(kcldtop)    -  hlf*qc(kcldtop)   )
         tmp1 = dqs(kcldtop)*tmp1/MAPL_CP

         tmp2 = ( qv(i,j,kcldtop-1)   +  qc(kcldtop-1) ) - &
                ( qv(i,j,kcldtop)     +  qc(kcldtop)   )  

         chis = -qc(kcldtop)*( 1 + hlf * dqs(kcldtop) / MAPL_CP )

         if ( ( tmp2 - tmp1 ) >= 0.0 ) then
            chis = 0.
         else
            chis = chis / ( tmp2 - tmp1 ) 
         endif

         if ( chis .gt. 1.0 ) chis=1.0

         slmixture = ( 1.0-chis )* ( slv(kcldtop)    -  hlf*qc(kcldtop)   )   &
                   +       chis  * ( slv(kcldtop-1)  -  hlf*qc(kcldtop-1) )


!-----------------------------------------------------------
! compute temperature of parcel at cloud top, svpcp.
         svpcp = slmixture /MAPL_CP

         buoypert   = ( slmixture - slv(kcldtop) )/MAPL_CP

!-----------------------------------------------------------
! calculate my best guess at the LCs' D parameter attributed 
! to Siems et al.
         stab       =    slv(kcldtop-1) - slv(kcldtop)
         if (stab .eq. 0.) then 
            dsiems  =  ( slv(kcldtop) - slmixture ) ! / 1.0  ! arbitrary, needs to be re-thought 
         else
            dsiems  =  ( slv(kcldtop) - slmixture ) / stab
         endif
         dsiems     =  min( dsiems, 10. )
         dsiems     =  max( dsiems,  0. )
         radf = maxradf
         zradtop = zhalf(i,j,kcldtop)

!-----------------------------------------------------------
! find depth of radiatively driven convection 

!-----------------------------------------------------------
! Expose radperturb and other funny business outside of radml_depth
         radperturb   = min( maxradf/100. , 0.3 ) ! dim. argument based on 100m deep cloud over 1000s
         do_jump_exit = .true.
         critjump     = 0.3
         svpcp = svpcp - radperturb

         slvcp = slv/MAPL_CP

         if (kcldtop .lt. ibot) then 
            call radml_depth(               &
                  i,j,icol,jcol,            &
                  nlev,kcldtop,ibot,        &
                  svpcp,zradtop,            &
                  critjump, do_jump_exit,   &
                  slvcp,                    &
                  zfull,                    &
                  zhalf,zradbase,           &
                  zradml                    )      
         else
            zradbase(i,j) = 0.0
            zradml(i,j)   = zradtop  
         end if

         zcloud(i,j) = zhalf(i,j,kcldtop) - zhalf(i,j,kcldbot)

         if (zradml(i,j) .le. 0.0 ) then 
            if(associated(radrcode_diag)) radrcode_diag(i,j)=3.
            go to 55   ! break out here if zradml<=0.0
         endif

!-----------------------------------------------------------
! compute radiation driven scale
!
! Vrad**3 = g*zradml*radf/density/slv

         vrad3 = MAPL_GRAV*zradml(i,j)*maxradf/density(kcldtop)/slv(kcldtop)   


!-----------------------------------------------------------
! compute entrainment rate
!

!-----------------------------------------------------------
! tmp1 here should be the buoyancy jump at cloud top
! SAK has it w/ resp to parcel property - svpcp. Im not 
! sure about that.
         tmp1 = MAPL_GRAV*max(0.1,((slv(kcldtop-1)/MAPL_CP)-svpcp))/(slv(kcldtop)/MAPL_CP)

!-----------------------------------------------------------
! Straightforward buoyancy jump across cloud top
         tmp1 = MAPL_GRAV*max( 0.1, ( slv(kcldtop-1)-slv(kcldtop) )/MAPL_CP ) &
                             / ( slv(kcldtop) /MAPL_CP )

!-----------------------------------------------------------
! compute buoy rev driven scale
         vbr3  = ( max( tmp1*zcloud(i,j), 0.)**3 )
         vbr3  = Abuoy*(chis**2)*max(dsiems,0.)*SQRT( vbr3 ) 

!----------------------------------------
! adjust velocity scales to prevent jumps 
! near zradtop=zcldtopmax
         if ( zradtop .gt. zcldtopmax-500. ) then 
            vrad3 = vrad3*(zcldtopmax - zradtop)/500.  
            vbr3  = vbr3 *(zcldtopmax - zradtop)/500.  
         endif
         vrad3=max( vrad3, 0. ) ! these really should not be needed
         vbr3 =max( vbr3,  0. )
!-----------------------------------------



         vrad = vrad3 ** (1./3.)    
         vbrv = vbr3  ** (1./3.)


         tmp2 = (  vrad**2 + vbrv**2  ) / zradml(i,j)
         wentr_rad = min(wentrmax,beta_rad*(vrad3+vbr3)/zradml(i,j)/(tmp1+tmp2))

         wentr_brv = beta_rad*vbr3/zradml(i,j)/(tmp1+tmp2)


!----------------------------------------
! fudgey adjustment of entrainment to reduce it
! for shallow boundary layers, and increase for 
! deep ones

!!AMM107
         if ( zradtop .lt. 500. ) then
            wentr_rad = 0.00
         endif
         if (( zradtop .gt. 500.) .and. (zradtop .le. 800. )) then
            wentr_rad = wentr_rad * ( zradtop-500.) / 300.
         endif

         if ( zradtop .lt. 2400. ) then 
            wentr_rad = wentr_rad * ( zradtop / 800. )
         else
            wentr_rad = 3.*wentr_rad
         endif
!-----------------------------------------

         k_entr_tmp = min ( akmax, wentr_rad*(zfull(i,j,kcldtop-1)-zfull(i,j,kcldtop)) )

         radfq(kcldtop)        = 1.
         radpbl                = 1.
         k_rad_col(kcldtop)    = k_entr_tmp
         k_t_entr(i,j,kcldtop) = k_t_entr(i,j,kcldtop) + k_entr_tmp
         k_m_entr(i,j,kcldtop) = k_m_entr(i,j,kcldtop) + k_entr_tmp


         if(associated(del_buoy_diag))  del_buoy_diag(i,j)    = tmp1
         if(associated(vrad_diag))      vrad_diag(i,j)        = vrad
         if(associated(kentrad_diag))   kentrad_diag(i,j)     = k_entr_tmp

         if(associated(chis_diag))      chis_diag(i,j)        = chis
         if(associated(vbrv_diag))      vbrv_diag(i,j)        = vbrv
         if(associated(dsiems_diag))    dsiems_diag(i,j)      = dsiems
         if(associated(wentr_brv_diag)) wentr_brv_diag(i,j)   = wentr_brv
         if(associated(wentr_rad_diag)) wentr_rad_diag(i,j)   = wentr_rad

         if(associated(zcldtop_diag))   zcldtop_diag(i,j)     = zhalf(i,j,kcldtop) 
         if(associated(slmixture_diag)) slmixture_diag(i,j)   = buoypert
         if(associated(delsinv_diag))   delsinv_diag(i,j)     = ( slv(kcldtop-1) - slv(kcldtop) )/MAPL_CP
         if(associated(cldradf_diag))   cldradf_diag(i,j)     = radf



!-----------------------------------------------------------
! handle case of radiatively driven top being the same top
! as surface driven top

         if (ipbl .eq. kcldtop .and. ipbl .gt. 0) then

            tmp2 = ((vbr3+vrad3+vsurf3+vshear3)**(2./3.)) / zradml(i,j)

            wentr_rad = min( wentrmax,  max(0., &
                  ((beta_surf *(vsurf3 + vshear3)+beta_rad*(vrad3+vbr3) )/ &
                  zradml(i,j))/(tmp1+tmp2) ) )

            wentr_pbl = wentr_rad

            k_entr_tmp = min ( akmax, wentr_rad*(zfull(i,j,kcldtop-1)-zfull(i,j,kcldtop)) )

            pblfq(ipbl)           = 1.
            radfq(kcldtop)        = 1.
            radpbl                = 1.
            k_rad_col(kcldtop)    = k_entr_tmp
            k_t_troen(ipbl)       = k_entr_tmp
            k_m_troen(ipbl)       = k_entr_tmp
            k_t_entr(i,j,kcldtop) = k_entr_tmp
            k_m_entr(i,j,kcldtop) = k_entr_tmp

         end if

!-----------------------------------------------------------
! if there are any interior layers to calculate diffusivity

         if ( kcldtop .lt. ibot ) then   

            do k = kcldtop+1,ibot

               ztmp = max(0.,(zhalf(i,j,k)-zradbase(i,j))/zradml(i,j) )

               if (ztmp.gt.0.) then

                  radfq(k) = 1.
                  k_entr_tmp = khradfac*MAPL_KARMAN*( vrad+vbrv )*ztmp* &
                        zradml(i,j)*ztmp*((1.-ztmp)**0.5)
                  k_entr_tmp = min ( k_entr_tmp, akmax )
                  k_rad_col(k) = k_entr_tmp
                  k_t_entr(i,j,k) = k_t_entr(i,j,k) + k_entr_tmp
                  k_m_entr(i,j,k) = k_m_entr(i,j,k) + k_entr_tmp*prandtlrad

               end if
            enddo

         end if

!-----------------------------------------------------------
! handle special case of zradbase < zsml
!
! in this case there should be no entrainment from the 
! surface.

         if (zradbase(i,j) .lt. zsml(i,j) .and. convpbl .eq. 1. .and. ipbl .gt. kcldtop) then
            wentr_pbl          = 0.
            pblfq(ipbl)        = 0.
            k_t_entr(i,j,ipbl) = k_t_entr(i,j,ipbl) - k_t_troen(ipbl)
            k_m_entr(i,j,ipbl) = k_m_entr(i,j,ipbl) - k_m_troen(ipbl)          
            k_t_troen(ipbl)    = 0.
            k_m_troen(ipbl)    = 0. 
         end if

55       continue


!-----------------------------------------------------------
!
! Modify diffusivity coefficients using MAX( A , B )        
         do k = 2, ibot     
            diff_t(i,j,k)   = max( k_t_entr(i,j,k), diff_t(i,j,k) )
            diff_m(i,j,k)   = max( k_m_entr(i,j,k), diff_m(i,j,k) )
            k_t_entr(i,j,k) = max( k_t_entr(i,j,k),           0.0 )
            k_m_entr(i,j,k) = max( k_m_entr(i,j,k),           0.0 )
         enddo

         do k = 1, nlev
            k_sfc(i,j,k) = k_t_troen(k)
            k_rad(i,j,k) = k_rad_col(k)
         end do

#ifndef _CUDA
      end do J_LOOP
      end do I_LOOP
#else
      end if J_LOOP
      end if I_LOOP
#endif



!-----------------------------------------------------------------------
! 
!      subroutine end
!

   end subroutine entrain


!======================================================================= 
!
!  Subroutine to calculate pbl depth
!
#ifdef _CUDA
   attributes(device) &
#endif
   subroutine mpbl_depth(i,j,icol,jcol,nlev,tpfac, entrate, pceff, vscale, pertopt, t, q, u, v, z, p, b_star, u_star , evap, sh, ipbl, ztop )

!
!  -----
!  INPUT
!  -----
!
!  i             column number
!  j             column number
!  icol          total number of columns
!  jcol          total number of columns
!  nlev          number of levels
!  t             temperature (K)
!  q             specific humidity (g/g)
!  u             zonal wind (m/s)
!  v             meridional wind (m/s)
!  b_star        buoyancy scale (m s-2)
!  u_star        surface velocity scale (m/s)
!       
!  ------
!  OUTPUT
!  ------
!
!  ipbl          half level containing pbl height
!  ztop          pbl height (m)

      integer, intent(in   )                            :: i, j, nlev, icol, jcol
      real,    intent(in   ), dimension(icol,jcol,nlev) :: t, z, q, p, u, v
      real,    intent(in   ), dimension(icol,jcol)      :: b_star, u_star, evap, sh
      real,    intent(in   )                            :: tpfac, entrate, pceff, vscale, pertopt
      integer, intent(  out)                            :: ipbl
      real,    intent(  out),dimension(icol,jcol)       :: ztop


      real     :: tep,z1,z2,t1,t2,qp,pp,qsp,dqp,dqsp,u1,v1,u2,v2,du
      real     :: entfr,entrate_x,lts,zrho,buoyflx,delzg,wstar
      integer  :: k


      !real, dimension(nlev) :: qst ! Not used in this code?


!calculate surface parcel properties

    if (pertopt /= 0) then
      zrho = p(i,j,nlev)/(287.04*(t(i,j,nlev)*(1.+0.608*q(i,j,nlev))))

      buoyflx = (sh(i,j)/MAPL_CP+0.608*t(i,j,nlev)*evap(i,j))/zrho ! K m s-1                                                                                                  
      delzg = (50.0)*MAPL_GRAV   ! assume 50m surface scale                                                                                                               
      wstar = max(0.,0.001+0.41*buoyflx*delzg/t(i,j,nlev)) ! m3 s-3      

      if (wstar > 0.001) then
        wstar = 1.0*wstar**.3333
!        print *,'sh=',sh(i,j),'evap=',evap(i,j),'wstar=',wstar
        tep  = t(i,j,nlev) + 0.4 + 2.*sh(i,j)/(zrho*wstar*MAPL_CP)
        qp   = q(i,j,nlev) + 2.*evap(i,j)/(zrho*wstar)
!        print *,'tpert=',2.*sh(i,j)/(zrho*wstar*MAPL_CP)
      else

      end if
    else   ! tpfac scales up bstar by inv. ratio of
           ! heat-bubble area to stagnant area
      if (nlev.eq.72) then
        tep  = (t(i,j,nlev) + 0.4) * (1.+ tpfac * b_star(i,j)/MAPL_GRAV)
      else
        tep  = (t(i,j,nlev) + 0.4) * (1.+ min(0.01,tpfac * b_star(i,j)/MAPL_GRAV))
      end if
      qp   = q(i,j,nlev)
    end if

!--------------------------------------------
! wind dependence of plume character. 
! 
!    actual_entrainment_rate_at_z  ~ entrate * [ 1.0 +  |U(z)-U(0)| / vscale ]
! 
! entrate:  tunable param from rc file
! vscale:   tunable param hardwired here.

! vscale is vscale_surf=0.25/100.0 parameter now passed through argument list

!search for level where this is exceeded              

      lts =  0.0
!  LTS using TH at 3km abve surface
      if (nlev.ne.72) then
         do k = nlev-1,2,-1
            if (z(i,j,k).gt.3000.0) then
              lts = t(i,j,k-1)*(1e5/p(i,j,k))**0.286
              exit
            end if
         end do
         lts = lts - t(i,j,nlev-1)*(1e5/p(i,j,nlev-1))**0.286
      end if

      t1   = t(i,j,nlev)
      v1   = v(i,j,nlev)
      u1   = u(i,j,nlev)
      z1   = z(i,j,nlev)
      ztop(i,j) = z1
      do k = nlev-1 , 2, -1
         z2 = z(i,j,k)
         t2 = t(i,j,k)
         u2 = u(i,j,k)
         v2 = v(i,j,k)
         pp = p(i,j,k)

!!Old Shear     du = sqrt ( ( u(i,j,k) - u1 )**2 + ( v(i,j,k) - v1 )**2 )
         du = sqrt ( ( u2 - u1 )**2 + ( v2 - v1 )**2 ) / (z2-z1)
         du = min(du,1.0e-8)

         entrate_x = entrate * ( 1.0 + du / vscale )

         entfr = min( entrate_x*(z2-z1), 0.99 )

         qp    = qp  + entfr*(q(i,j,k)-qp)

! dry adiabatic ascent through one layer.
! Static energy conserved. 
         tep   = tep - MAPL_GRAV*( z2-z1 )/MAPL_CP

! Environmental air entrained
         tep   = tep + entfr*(t(i,j,k)-tep)

         dqsp  = dqsat(tep , pp , qsat=qsp,  pascals=.true. )

         dqp   = max( qp - qsp, 0. )/(1.+(MAPL_ALHL/MAPL_CP)*dqsp )
         qp    = qp - dqp
         if (lts .eq. 0.0) then
           tep   = tep  + pceff * MAPL_ALHL * dqp/MAPL_CP  ! "Precipitation efficiency" basically means fraction
! of condensation heating that gets applied to parcel
         else
           tep = tep + (pceff + 0.5*(1.-pceff)*(1.+TANH(lts-18.)))*MAPL_ALHL * dqp/MAPL_CP      
!                           Set pceff to 1 where LTS is high
!           tep   = tep  + pceff * MAPL_ALHL * dqp/MAPL_CP  ! "Precipitation efficiency" basically means fraction
!                                                           ! of condensation heating that gets applied to parcel
         endif

! If parcel temperature (tep) colder than env (t2)
! OR if entrainment too big, declare this the PBL top
         if ( (t2 .ge. tep) .or. ( entfr .ge. 0.9899 ) ) then
!!a0082  Make the parcel a little less buoyant (liquid water loading)- 
!!a0082     translates here into subtracting from tep
!!a0082   if ( (t2 .ge. (tep-0.2) ) .or. ( entfr .ge. 0.9899 ) ) then
            ztop(i,j) = 0.5*(z2+z1)
            ipbl = k+1
            exit
         end if

         z1 = z2
         t1 = t2
         u1 = u2
         v1 = v2
      enddo


      return

   end subroutine mpbl_depth
!=======================================================================

!======================================================================= 
!
!  Subroutine to calculate bottom and depth of radiatively driven mixed
!  layer
!
!

#ifdef _CUDA
   attributes(device) &
#endif
   subroutine radml_depth(i, j, icol, jcol, nlev, toplev, botlev, &
         svp, zt, critjump, do_jump_exit, t, zf, zh, zb, zml)

!
!  -----
!  INPUT
!  -----
!
!  i        column number
!  j        column number
!  icol     total number of columns
!  jcol     total number of columns
!  toplev   top level of calculation
!  botlev   bottom level of calculation
!  nlev     number of levels
!  svp      cloud top liquid water virtual static energy divided by cp (K)
!  zt       top of radiatively driven layer (m)
!  t        liquid water virtual static energy divided by cp (K)
!  zf       full level height above ground (m)
!  zh       half level height above ground (m)
!       
!  ------
!  OUTPUT
!  ------
!
!  zb      base height of radiatively driven mixed layer (m)
!  zml     depth of radiatively driven mixed layer (m)


      integer, intent(in   )                              :: i, j, toplev, botlev, icol, jcol, nlev
      real,    intent(in   )                              :: svp, zt, critjump
      real,    intent(in   ), dimension(nlev)             :: t
      real,    intent(in   ), dimension(icol,jcol,nlev)   :: zf
      real,    intent(in   ), dimension(icol,jcol,nlev+1) :: zh
      real,    intent(  out), dimension(icol,jcol)        :: zb, zml
      logical, intent(in   )                              :: do_jump_exit

      real    :: svpar,h1,h2,t1,t2,entrate,entfr
      integer :: k


      !initialize zml
      zml(i,j) = 0.

      !calculate cloud top parcel properties
      svpar   = svp
      h1      = zf(i,j,toplev)
      t1      = t(toplev)
      if (nlev.eq.72) then
        entrate = 0.2/200.
      else
        entrate = 1.0/1000.
      endif

      !search for level where parcel is warmer than env             

      ! first cut out if parcel is already warmer than
      ! cloudtop. 
      if (t1.lt.svpar) then
         zb(i,j)  = h1
         zml(i,j) = 0.00
         return
      endif

      ! If above is false keep looking
      do k = 1,botlev
         if (k > toplev) then
            h2 = zf(i,j,k)
            t2 = t(k)

            if (t2.lt.svpar) then
               if ( abs(t1 - t2 ) .gt. 0.2 ) then 
                  zb(i,j) = h2 + (h1 - h2)*(svpar - t2)/(t1 - t2)
                  zb(i,j) = MAX( zb(i,j) , 0. )  ! final protection against interp problems
               else
                  zb(i,j) = h2
               endif
               zml(i,j) = zt - zb(i,j)
               return
            end if

            if (do_jump_exit .and. (t1-t2) .gt. critjump .and. k .gt. toplev+1) then
               zb(i,j)  = zh(i,j,k)
               zml(i,j) = zt - zb(i,j)
               return
            end if

            entfr = min( entrate*(h1-h2), 1.0 )
            svpar = svpar + entfr*(t2-svpar)

            h1 = h2
            t1 = t2
         end if
      enddo

      zb(i,j)  = 0.
      zml(i,j) = zt

      return
   end subroutine radml_depth

!=======================================================================
!========================================================================  
!       Subroutine to return the vertical K-profile of diffusion 
!       coefficients for the surface driven convective mixed layer.
!       This code returns to form used Lock et al.. Does not match
!       to surface layer.  
!    
!   call diffusivity_pbl2(i, j, lm, h, k_ent, vsurf, zm, k_m, k_t)
!		
!      i:      Column number
!      j:      Column number
!      icol:   Total number of columns
!      jcol:   Total number of columns
!      lm:     Number of levels (one less than total in half-level zm)
!      h:      Depth of surface driven mixed layer (m) 
!      k_ent:  PBL top entrainment diffusivity (m+2 s-1)
!      vsurf:  PBL top entrainment velocity scale (m s-1)
!      zm:     Half level heights relative to the ground (m),        DIM[1:lm+1]
!      k_m:    Momentum diffusion coefficient (m+2 s-1),             DIM[1:lm] (edges)
!      k_t:    Heat and tracer diffusion coefficient (m+2 s-1)       DIM[1:lm] (edges)

#ifdef _CUDA
   attributes(device) &
#endif
   subroutine diffusivity_pbl2(i, j, icol, jcol, lm, h, kfac_l, kfac_o, k_ent, vsurf, frland, zm, k_m, k_t)

      integer, intent(in   )                              :: i, j, lm, icol, jcol
      real,    intent(in   ), dimension(icol,jcol)        :: h, frland
      real,    intent(in   )                              :: k_ent, vsurf, kfac_l, kfac_o
      real,    intent(in   ), dimension(icol,jcol,1:lm+1) :: zm
      real,    intent(  out), dimension(lm)               :: k_m, k_t

      real    :: EE, hin, kfacx 

      integer :: k

      kfacx = kfac_l*frland(i,j) + kfac_o*(1.0-frland(i,j))

      hin = 0.0 ! 200.  ! "Organization" scale for plume (m).

      do k = 1, lm
         k_m(k) = 0.0
         k_t(k) = 0.0
      end do

!! factor = (zm(i,j,k)/hinner)* (1.0 -(zm(i,j,k)-hinner)/(h(i,j)-hinner))**2

      if ( vsurf*h(i,j) .gt. 0. ) then
         EE  = 1.0 - sqrt( k_ent / ( kfacx * MAPL_KARMAN * vsurf * h(i,j) ) )
         EE  = max( EE , 0.7 )  ! If EE is too small, then punt, as LCs
         do k = 1, lm
            if ( ( zm(i,j,k) .le. h(i,j) ) .and.  ( zm(i,j,k) .gt. hin )  ) then
               k_t(k) = kfacx * MAPL_KARMAN * vsurf * ( zm(i,j,k)-hin ) * ( 1. - EE*( (zm(i,j,k)-hin)/(h(i,j)-hin) ))**2
            end if
         end do
      endif

      do k = 1, lm
         k_m(k) = k_t(k)
      end do

      return
   end subroutine diffusivity_pbl2

!
!======================================================================= 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _CUDA 
   attributes(device) function QSAT(TL,PL,PASCALS)

     real,              intent(IN) :: TL, PL
     logical, optional, intent(IN) :: PASCALS
     real    :: QSAT

     real    :: URAMP, DD, QQ, TI, DQ, PP
     integer :: IT

     URAMP = TMIX

     if (present(PASCALS)) then 
        if (PASCALS) then
           PP = PL
        else
           PP = PL*100.
        end if
     else
        PP = PL*100.
     end if

     TI = TL - ZEROC

     if    (TI <= URAMP) then
        QSAT  =  QSATICE0(TL,PP,DQ) 
     elseif(TI >= 0.0  ) then
        QSAT  =  QSATLQU0(TL,PP,DQ)
     else
        QSAT  =  QSATICE0(TL,PP,DQ)
        QQ    =  QSATLQU0(TL,PP,DQ)
        TI    =  TI/URAMP
        QSAT  =  TI*(QSAT - QQ) +  QQ
     end if

   end function QSAT

   attributes(device) function DQSAT(TL,PL,QSAT,PASCALS)

      real,              intent(IN) :: TL, PL
      real,              intent(OUT):: QSAT
      logical, optional, intent(IN ):: PASCALS
      real    :: DQSAT
      real    :: URAMP, TT, WW, DD, DQQ, QQ, TI, DQI, QI, PP, DQ
      integer :: IT

      URAMP = TMIX

      if (present(PASCALS)) then 
         if (PASCALS) then
            PP = PL
         else
            PP = PL*100.
         end if
      else
         PP = PL*100.
      end if

      TI = TL - ZEROC

      if    (TI <= URAMP) then
         QQ  = QSATICE0(TL,PP,DQ)
         QSAT  = QQ
         DQSAT = DQ
      elseif(TI >= 0.0  ) then
         QQ  = QSATLQU0(TL,PP,DQ)
         QSAT  = QQ
         DQSAT = DQ
      else
         QQ  = QSATLQU0(TL,PP,DQQ)
         QI  = QSATICE0(TL,PP,DQI)
         TI  = TI/URAMP
         DQSAT = TI*(DQI - DQQ) + DQQ
         QSAT  = TI*(QI - QQ) +  QQ
      end if

   end function DQSAT

   attributes(device) function QSATLQU0(TL,PL,DQ) result(QS)

      real, intent(IN) :: TL
      real, intent(IN) :: PL
      real, intent(OUT):: DQ
      real    :: QS
   
      real    :: TI,W
      real    :: DD
      real    :: TT
      real    :: DDQ
      integer :: IT

      integer, parameter :: TYPE = 1

#define TX TL 
#define PX PL 
#define EX QS
#define DX DQ 


   if    (TX<TMINLQU) then
      TI = TMINLQU
   elseif(TX>TMAXTBL) then
      TI = TMAXTBL 
   else
      TI = TX  
   end if

#include "esatlqu.code"

   if    (TX<TMINLQU) then
      DDQ = 0.0
   elseif(TX>TMAXTBL) then
      DDQ = 0.0
   else
      if(PX>EX) then
         DD = EX 
         TI = TX + DELTA_T
#include "esatlqu.code"
         DDQ = EX-DD
         EX  = DD
      endif
   end if

   if(PX > EX) then
      DD = ESFAC/(PX - (1.0-ESFAC)*EX)
      EX = EX*DD
      DX = DDQ*ERFAC*PX*DD*DD
   else
      EX = MAX_MIXING_RATIO
      DX = 0.0
   end if

#undef  DX
#undef  TX
#undef  EX
#undef  PX

      return
   end function QSATLQU0

   attributes(device) function QSATICE0(TL,PL,DQ) result(QS)

      real, intent(IN) :: TL
      real, intent(IN) :: PL
      real, intent(OUT):: DQ
      real    :: QS

      real    :: TI,W
      real    :: DD
      real    :: TT
      real    :: DDQ
      integer :: IT

      integer, parameter :: TYPE = 1

#define TX TL
#define PX PL
#define EX QS
#define DX DQ


   if    (TX<TMINICE) then
      TI = TMINICE
   elseif(TX>ZEROC  ) then
      TI = ZEROC
   else
      TI = TX
   end if

#include "esatice.code"

   if    (TX<TMINICE) then
      DDQ = 0.0
   elseif(TX>ZEROC  ) then
      DDQ = 0.0
   else
      if(PX>EX) then
         DD = EX
         TI = TX + DELTA_T
#include "esatice.code"
         DDQ = EX-DD
         EX  = DD
      endif
   end if

   if(PX > EX) then
      DD = ESFAC/(PX - (1.0-ESFAC)*EX)
      EX = EX*DD
      DX = DDQ*ERFAC*PX*DD*DD
   else
      EX = MAX_MIXING_RATIO
      DX = 0.0
   end if

#undef  DX
#undef  TX
#undef  EX
#undef  PX

         return
   end function QSATICE0

#endif

!*********************************************************************

end module LockEntrain
