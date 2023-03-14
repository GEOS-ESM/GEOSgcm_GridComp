module uwshcu

!#define UWDIAG 1

   use GEOS_Mod, only: write_parallel
   use MAPL,     only: MAPL_UNDEF

   use GEOS_UtilsMod, only: GEOS_QSAT, GEOS_DQSAT
   use GEOSmoist_Process_Library
   use MAPL_ConstantsMod, only: MAPL_TICE , MAPL_CP   , &
                                MAPL_GRAV , MAPL_ALHS , &
                                MAPL_ALHL , MAPL_ALHF , &
                                MAPL_RGAS , MAPL_H2OMW, &
                                MAPL_AIRMW, MAPL_RVAP , &
                                MAPL_PI, r8 => MAPL_R8

   implicit none

  type SHLWPARAM_TYPE
     integer  :: niter_xc           ! Number xc iterations
     integer  :: iter_cin           ! Number iterations for implicit CIN
     integer  :: use_CINcin         ! if true, calc CIN thru L..
     integer  :: use_self_detrain   ! 
     integer  :: use_momenflx       ! Perform momentum transport
     integer  :: use_cumpenent      ! Cumulative penetrative entrainment
     integer  :: scverbose          ! activate print statements
     integer  :: windsrcavg         ! Source air uses PBL mean momentum
     real     :: rpen               ! Penentrative entrainment factor
     real     :: rle
     real     :: rkm                ! Factor controlling lateral mixing rate
     real     :: mixscale           ! Controls vertical structure of mixing
     real     :: detrhgt            ! Mixing rate increases above this height
     real     :: rkfre              ! Vertical velocity fraction of tke
     real     :: rmaxfrac           ! Maximum core updraft fraction
     real     :: mumin1             ! 
     real     :: rbuoy              ! Non-hydro pressure effect on updraft
     real     :: rdrag              ! Drag coefficient
     real     :: epsvarw            ! Variance of PBL w by mesoscale
     real     :: PGFc               ! Pressure gradient force
     real     :: criqc              ! Updraft maximum condensate 
     real     :: frc_rasn           ! Precip fraction of expelled condensate
     real     :: kevp               ! Evaporative efficiency
     real     :: rdrop              ! liquid drop radius
     real     :: thlsrc_fac         ! Scaling factor for thlsrc perturbation
     real     :: qtsrc_fac          ! Scaling factor for qtsrc perturbation
     real     :: qtsrchgt           ! Interpolation height for total water source
     integer  :: cridist_opt
  endtype SHLWPARAM_TYPE
  type   (SHLWPARAM_TYPE) :: SHLWPARAMS

   private
   public compute_uwshcu_inv, SHLWPARAMS

   real, parameter :: xlv   = MAPL_ALHL             ! Latent heat of vaporization
   real, parameter :: xlf   = MAPL_ALHF             ! Latent heat of fusion
   real, parameter :: xls   = MAPL_ALHS             ! Latent heat sublimation
   real, parameter :: cp    = MAPL_CP               ! Specific heat of dry air
   real, parameter :: zvir  = 0.609                 ! r_H2O/r_air-1
   real, parameter :: r     = MAPL_RGAS             ! Gas constant for dry air
   real, parameter :: g     = MAPL_GRAV             ! Acceleration due to gravity
   real, parameter :: ep2   = MAPL_H2OMW/MAPL_AIRMW ! Ratio of molc wgt H2O to dry air
   real, parameter :: p00   = 1e5                   ! Reference pressure
   real, parameter :: rovcp = MAPL_RGAS/MAPL_CP     ! Gas constant over specific heat

   real, parameter :: mintracer = tiny(1.)
contains

   real function exnerfn(pressure)
      real, intent(in)              :: pressure   ! in Pa
      exnerfn = (pressure/p00)**rovcp
      return
   end function exnerfn

   REAL FUNCTION fract_liq_f(temp2) ! temp2 in Kelvin, fraction between 0 and 1.
   implicit none
   real,intent(in) :: temp2 ! K
   real :: temp,ptc
   real, parameter :: max_temp = 46. !Celsius
!   SELECT CASE(FRAC_MODIS)
!   CASE (1)
   temp = temp2-273.16 !Celsius
   temp = min(max_temp,max(-max_temp,temp))
   ptc  = 7.6725 + 1.0118*temp + 0.1422*temp**2 + &
             0.0106*temp**3 + 3.39e-4 * temp**4 + &
             3.95e-6 * temp**5
   fract_liq_f = 1./(1.+exp(-ptc))
!   CASE DEFAULT
!   fract_liq_f = min(1., (max(0.,(temp2-t_ice))/(t_0-t_ice))**2)
!   END SELECT
   END FUNCTION


   subroutine compute_uwshcu_inv(idim, k0,        dt,pmid0_inv,     & ! INPUT
         zmid0_inv, exnmid0_inv, pifc0_inv, zifc0_inv, exnifc0_inv, &
         dp0_inv, u0_inv, v0_inv, qv0_inv, ql0_inv, qi0_inv,        &
         th0_inv, tke_inv, kpbl_inv, shfx,evap, cnvtr, frland,      & 
         cush,                                                      & ! INOUT
         umf_inv, dcm_inv, qvten_inv, qlten_inv, qiten_inv, thten_inv, & ! OUTPUT
         uten_inv, vten_inv, qrten_inv, qsten_inv, cufrc_inv,       &
         fer_inv, fdr_inv, qldet_inv, qidet_inv, qlsub_inv,         &
         qisub_inv, ndrop_inv, nice_inv, tpert_out, qpert_out,      & 
         qtflx_inv, slflx_inv, uflx_inv, vflx_inv,                  &
#ifdef UWDIAG
         qcu_inv, qlu_inv, qiu_inv, cbmf, qc_inv,                   & ! DIAGNOSTIC ONLY
         cnt_inv, cnb_inv, cin, plcl, plfc, pinv, prel, pbup,       &
         wlcl, qtsrc, thlsrc, thvlsrc, tkeavg, cldtop, wu_inv,      &
         qtu_inv, thlu_inv, thvu_inv, uu_inv, vu_inv, xc_inv,       &
#endif
         dotransport)

      implicit none


      integer, intent(in)   :: idim                     ! Number of columns
      integer, intent(in)   :: k0                       ! Number of levels  
      integer, intent(in)   :: dotransport              ! Transport tracers [1 true]
      real   , intent(in)   :: dt                       ! moist heartbeat [s]

      real,   intent(in)    :: pifc0_inv(idim,k0+1)     !  Environmental pressure at the interfaces [ Pa ]
      real,   intent(in)    :: zifc0_inv(idim,k0+1)     !  Environmental height at the interfaces   [ m ]
      real,   intent(in)    :: exnifc0_inv(idim,k0+1)   !  Exner function at the interfaces
      real,   intent(in)    :: pmid0_inv(idim,k0)       !  Environmental pressure at the layer mid-point [ Pa ]
      real,   intent(in)    :: zmid0_inv(idim,k0)       !  Environmental height at the layer mid-point [ m ]
      real,   intent(in)    :: exnmid0_inv(idim,k0)     !  Exner function at the layer mid-point
      real,   intent(in)    :: dp0_inv(idim,k0)         !  Environmental layer pressure thickness [ Pa ] > 0.
      real,   intent(in)    :: u0_inv(idim,k0)          !  Environmental zonal wind [ m/s ]
      real,   intent(in)    :: v0_inv(idim,k0)          !  Environmental meridional wind [ m/s ]
      real,   intent(in)    :: qv0_inv(idim,k0)         !  Environmental water vapor specific humidity [ kg/kg ]
      real,   intent(in)    :: ql0_inv(idim,k0)         !  Environmental liquid water specific humidity [ kg/kg ]
      real,   intent(in)    :: qi0_inv(idim,k0)         !  Environmental ice specific humidity [ kg/kg ]
      real,   intent(in)    :: th0_inv(idim,k0)         !  Environmental temperature [ K ]
      real,   intent(in)    :: tke_inv(idim,k0+1)       !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
                                                        !  at the previous time step [ fraction ]
      real, intent(in)    :: kpbl_inv(idim)           !  Height of PBL [ m ]
      real, intent(in)    :: shfx(idim)               ! Surface sensible heat
      real, intent(in)    :: evap(idim)               ! Surface evaporation
      real, intent(in)    :: cnvtr(idim)              ! convective tracer
      real, intent(in)    :: frland(idim)             ! land fraction
      real, intent(inout) :: cush(idim)               !  Convective scale height [m]

      real, intent(out)   :: umf_inv(idim,k0+1)       !  Updraft mass flux at interfaces [kg/m2/s]
      real, intent(out)   :: dcm_inv(idim,k0)       !  Detrained cloudy air mass
      real, intent(out)   :: qvten_inv(idim,k0)       !  Tendency of water vapor specific humidity [ kg/kg/s ]
      real, intent(out)   :: qlten_inv(idim,k0)         !  Tendency of liquid water specific humidity [ kg/kg/s ]
      real, intent(out)   :: qiten_inv(idim,k0)         !  Tendency of ice specific humidity [ kg/kg/s ]
      real, intent(out)   :: thten_inv(idim,k0)       !  Tendency of potential temperature [ K/s ]
      real, intent(out)   :: uten_inv(idim,k0)        !  Tendency of zonal wind [ m/s2 ]
      real, intent(out)   :: vten_inv(idim,k0)        !  Tendency of meridional wind [ m/s2 ]
!      real, intent(out)   :: trten_inv(idim,k0,ncnst) !  Tendency of tracers [ #/s, kg/kg/s ]
      real, intent(out)   :: qrten_inv(idim,k0)       !  Tendency of rain water specific humidity [ kg/kg/s ]
      real, intent(out)   :: qsten_inv(idim,k0)       !  Tendency of snow specific humidity [ kg/kg/s ]
      real, intent(out)   :: cufrc_inv(idim,k0)       !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
      real, intent(out)   :: fer_inv(idim,k0)
      real, intent(out)   :: fdr_inv(idim,k0)
      real, intent(out)   :: qldet_inv(idim,k0)
      real, intent(out)   :: qidet_inv(idim,k0)
      real, intent(out)   :: qlsub_inv(idim,k0)
      real, intent(out)   :: qisub_inv(idim,k0)
      real, intent(out)   :: ndrop_inv(idim,k0)
      real, intent(out)   :: nice_inv(idim,k0)
      real, intent(out)   :: tpert_out(idim)
      real, intent(out)   :: qpert_out(idim)
      real, intent(out)   :: qtflx_inv(idim,k0+1)
      real, intent(out)   :: slflx_inv(idim,k0+1)
      real, intent(out)   :: uflx_inv(idim,k0+1)
      real, intent(out)   :: vflx_inv(idim,k0+1)

!!! Diagnostic only
#ifdef UWDIAG
      real, intent(out)   :: qcu_inv(idim,k0)         !  Liquid+ice specific humidity within cumulus updraft [ kg/kg ]
      real, intent(out)   :: qlu_inv(idim,k0)         !  Liquid water specific humidity within cumulus updraft [ kg/kg ]
      real, intent(out)   :: qiu_inv(idim,k0)         !  Ice specific humidity within cumulus updraft [ kg/kg ]
      real, intent(out)   :: qc_inv(idim,k0)          !  Tendency of cumulus condensate detrained into the environment [ kg/kg/s ]
      real, intent(out)   :: wu_inv(idim,k0+1)
      real, intent(out)   :: qtu_inv(idim,k0+1)
      real, intent(out)   :: thlu_inv(idim,k0+1)
      real, intent(out)   :: thvu_inv(idim,k0+1)
      real, intent(out)   :: uu_inv(idim,k0+1)
      real, intent(out)   :: vu_inv(idim,k0+1)
      real, intent(out)   :: xc_inv(idim,k0)
      real, intent(out)   :: cbmf(idim)                !  Cumulus base mass flux [ kg/m2/s ]
      real, intent(out)   :: cnt_inv(idim)             !  Cumulus top  interface index, cnt = kpen [ no ]
      real, intent(out)   :: cnb_inv(idim)             !  Cumulus base interface index, cnb = krel - 1 [ no ]
      real, intent(out)   :: cin(idim)
      real, intent(out)   :: plcl(idim)
      real, intent(out)   :: plfc(idim)
      real, intent(out)   :: pinv(idim)
      real, intent(out)   :: prel(idim)
      real, intent(out)   :: pbup(idim)
      real, intent(out)   :: wlcl(idim)
      real, intent(out)   :: qtsrc(idim)
      real, intent(out)   :: thlsrc(idim)
      real, intent(out)   :: thvlsrc(idim)
      real, intent(out)   :: tkeavg(idim)
      real, intent(out)   :: cldtop(idim)
#endif

  !----- Local variables -----
      real              :: pifc0(idim,0:k0)         !  Environmental pressure at the interfaces [ Pa ]
      real              :: zifc0(idim,0:k0)         !  Environmental height at the interfaces   [ m ]
      real              :: exnifc0(idim,0:k0)       !  Exner function on interfaces
      real              :: pmid0(idim,k0)           !  Environmental pressure at the layer mid-point [ Pa ]
      real              :: zmid0(idim,k0)           !  Environmental height at the layer mid-point [ m ]
      real              :: exnmid0(idim,k0)         !  Exner function on layer mid-point
      real              :: dp0(idim,k0)             !  Environmental layer pressure thickness [ Pa ] > 0.
      real              :: u0(idim,k0)              !  Environmental zonal wind [ m/s ]
      real              :: v0(idim,k0)              !  Environmental meridional wind [ m/s ]
      real              :: qv0(idim,k0)             !  Environmental water vapor specific humidity [ kg/kg ]
      real              :: ql0(idim,k0)             !  Environmental liquid water specific humidity [ kg/kg ]
      real              :: qi0(idim,k0)             !  Environmental ice specific humidity [ kg/kg ]
      real              :: th0(idim,k0)             !  Environmental temperature [ K ]
      real              :: tke(idim,0:k0)           !  Turbulent kinetic energy [ m2 s-2 ] 
      real, allocatable :: tr0(:,:,:)               !  Environmental tracers [ #, kg/kg ]
      real              :: umf(idim,0:k0)           !  Updraft mass flux at the interfaces [ kg/m2/s ]
      real              :: dcm(idim,k0)             !  Detrained cloudy air mass
      real              :: qvten(idim,k0)           !  Tendency of water vapor specific humidity [ kg/kg/s ]
      real              :: qlten(idim,k0)           !  Tendency of liquid water specific humidity [ kg/kg/s ]
      real              :: qiten(idim,k0)           !  tendency of ice specific humidity [ kg/kg/s ]
      real              :: sten(idim,k0)            !  Tendency of static energy [ J/kg/s ]
      real              :: uten(idim,k0)            !  Tendency of zonal wind [ m/s2 ]
      real              :: vten(idim,k0)            !  Tendency of meridional wind [ m/s2 ]
      real              :: qrten(idim,k0)           !  Tendency of rain water specific humidity [ kg/kg/s ]
      real              :: qsten(idim,k0)           !  Tendency of snow speficif humidity [ kg/kg/s ]
      real              :: cufrc(idim,k0)           !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
      real              :: fer(idim,k0)             !  Fractional entrainment rate
      real              :: fdr(idim,k0)             !  Fractional detrainment rate
      real              :: qldet(idim,k0)           !  Detrained liquid condensate
      real              :: qidet(idim,k0)           !  Detrained ice
      real              :: qlsub(idim,k0)           !  Liquid water tendency due to subsidence
      real              :: qisub(idim,k0)           !  Ice water tendency due to subsidence
      real              :: ndrop(idim,k0)
      real              :: nice(idim,k0)
      real              :: qtflx(idim,0:k0)
      real              :: slflx(idim,0:k0)
      real              :: uflx(idim,0:k0)
      real              :: vflx(idim,0:k0)
      real              :: cnvtrmax(idim)
      real              :: tmp2d(idim,k0)

!--------- Local, Diagnostic only ---------
#ifdef UWDIAG
      real              :: trten(idim,k0,ncnst)    !  Tendency of tracers [ #/s, kg/kg/s ]
      real              :: qcu(idim,k0)            !  Condensate water specific humidity within cumulus updraft
                                                   ! at the layer mid-point [ kg/kg ]
      real              :: qlu(idim,k0)            !  Liquid water specific humidity within cumulus updraft
                                                   ! at the layer mid-point [ kg/kg ]
      real              :: qiu(idim,k0)            !  Ice specific humidity within cumulus updraft
                                                   ! at the layer mid-point [ kg/kg ]
      real              :: qc(idim,k0)             !  Tendency of cumulus condensate detrained into the environment [ kg/kg/s ]
      real              :: cnt(idim)               !  Cumulus top  interface index, cnt = kpen [ no ]
      real              :: cnb(idim)               !  Cumulus base interface index, cnb = krel - 1 [ no ] 
      real              :: wu(idim,0:k0)
      real              :: qtu(idim,0:k0)
      real              :: thlu(idim,0:k0)
      real              :: thvu(idim,0:k0)
      real              :: uu(idim,0:k0)
      real              :: vu(idim,0:k0)
      real              :: xc(idim,k0)
      real              :: trten_inv(idim,k0,ncnst) !  Tendency of tracers [ #/s, kg/kg/s ]
#endif



!---------- Indices -----------
      integer           :: i                        !  Horizontal index for local fields [ no ] 
      integer           :: k                        !  Vertical index for local fields [ no ] 
      integer           :: k_inv                    !  Vertical index for incoming fields [ no ]
      integer           :: m                        !  Tracer index [ no ]
      integer           :: kpbl(idim)
      integer           :: ncnst, IM, JM

      ncnst = size(CNV_Tracers)
      IM = size(CNV_Tracers(1)%Q,1)
      JM = size(CNV_Tracers(1)%Q,2)
      allocate(tr0(idim,k0,ncnst))

      ! flip mid-level variables
      do k = 1, k0
         k_inv               = k0 + 1 - k
         pmid0(:idim,k)      = pmid0_inv(:idim,k_inv)
         u0(:idim,k)         = u0_inv(:idim,k_inv)
         v0(:idim,k)         = v0_inv(:idim,k_inv)
         zmid0(:idim,k)      = zmid0_inv(:idim,k_inv)
         exnmid0(:idim,k)    = exnmid0_inv(:idim,k_inv)
         dp0(:idim,k)        = dp0_inv(:idim,k_inv)
         qv0(:idim,k)        = qv0_inv(:idim,k_inv)
         ql0(:idim,k)        = ql0_inv(:idim,k_inv)
         qi0(:idim,k)        = qi0_inv(:idim,k_inv)
         th0(:idim,k)        = th0_inv(:idim,k_inv)
         do m = 1, ncnst
            tr0(:idim,k,m)   = reshape(CNV_Tracers(m)%Q(:,:,k_inv), (/idim/))
         enddo
      enddo

      ! flip interface variables
      tke(:,:) = 0.
      pifc0(:,:) = 0.
      zifc0(:,:) = 0.
      exnifc0(:,:) = 0.
      do k = 0, k0
         k_inv               = k0 - k + 1
         tke(:idim,k)        = tke_inv(:idim,k_inv)
         pifc0(:idim,k)      = pifc0_inv(:idim,k_inv)
         zifc0(:idim,k)      = zifc0_inv(:idim,k_inv)
         exnifc0(:idim,k)    = exnifc0_inv(:idim,k_inv)
      end do

      kpbl = int(kpbl_inv)

      do i = 1,idim
!        cnvtrmax(i) = min(300.,max(0.,maxval(cnvtr(i,:))))
        cnvtrmax(i) = min(1e-5,max(0.,cnvtr(i)))
        if (frland(i)>0.5) cnvtrmax(i) = 0.
        if (isnan(cnvtrmax(i))) cnvtrmax(i) = 0.
      end do

      call compute_uwshcu( idim,k0, dt, ncnst,pifc0, zifc0, &
           exnifc0, pmid0, zmid0, exnmid0, dp0, u0, v0,     &
           qv0, ql0, qi0, th0, tr0, kpbl, tke, cush, umf,   &
           dcm, qvten, qlten, qiten, sten, uten, vten,      &
           qrten, qsten, cufrc, fer, fdr, qldet, qidet,     & 
           qlsub, qisub, ndrop, nice,                       &
           shfx, evap, cnvtrmax, tpert_out, qpert_out,      &
           qtflx, slflx, uflx, vflx,                        &
#ifdef UWDIAG
           qcu, qlu, qiu, cbmf, qc, cnt, cnb,               & ! Diagnostic only
           cin, plcl, plfc, pinv, prel, pbup, wlcl, qtsrc,  &
           thlsrc, thvlsrc, tkeavg, cldtop, wu, qtu,        &
           thlu, thvu, uu, vu, xc, trten,                   & 
#endif
           dotransport )

      ! Reverse again

#ifdef UWDIAG
      cnt_inv(:idim) = k0 + 1 - cnt(:idim)
      cnb_inv(:idim) = k0 + 1 - cnb(:idim)
#endif

      do k = 0, k0
         k_inv                    = k0 + 1 - k
         umf_inv(:idim,k_inv)     = umf(:idim,k)       
         qtflx_inv(:idim,k_inv)   = qtflx(:idim,k)       
         slflx_inv(:idim,k_inv)   = slflx(:idim,k)       
         uflx_inv(:idim,k_inv)    = uflx(:idim,k)       
         vflx_inv(:idim,k_inv)    = vflx(:idim,k)       

#ifdef UWDIAG
         wu_inv(:idim,k_inv)      = wu(:idim,k)   ! Diagnostic only
         qtu_inv(:idim,k_inv)     = qtu(:idim,k)
         thlu_inv(:idim,k_inv)    = thlu(:idim,k)
         thvu_inv(:idim,k_inv)    = thvu(:idim,k)
         uu_inv(:idim,k_inv)      = uu(:idim,k)
         vu_inv(:idim,k_inv)      = vu(:idim,k)
#endif
      end do

      do k = 1, k0
         k_inv                    = k0 + 1 - k
         dcm_inv(:idim,k_inv)     = dcm(:idim,k)    
         qvten_inv(:idim,k_inv)   = qvten(:idim,k)   
         qlten_inv(:idim,k_inv)   = qlten(:idim,k)   
         qiten_inv(:idim,k_inv)   = qiten(:idim,k)   
         thten_inv(:idim,k_inv)   = sten(:idim,k) / (cp*exnmid0(:idim,k))
         uten_inv(:idim,k_inv)    = uten(:idim,k)    
         vten_inv(:idim,k_inv)    = vten(:idim,k)    
         qrten_inv(:idim,k_inv)   = qrten(:idim,k)   
         qsten_inv(:idim,k_inv)   = qsten(:idim,k)  
         cufrc_inv(:idim,k_inv)   = cufrc(:idim,k)
         fer_inv(:idim,k_inv)     = fer(:idim,k)
         fdr_inv(:idim,k_inv)     = fdr(:idim,k)
         qldet_inv(:idim,k_inv)   = qldet(:idim,k)
         qidet_inv(:idim,k_inv)   = qidet(:idim,k)
         qlsub_inv(:idim,k_inv)   = qlsub(:idim,k)
         qisub_inv(:idim,k_inv)   = qisub(:idim,k)
         ndrop_inv(:idim,k_inv)   = ndrop(:idim,k)
         nice_inv(:idim,k_inv)    = nice(:idim,k)
#ifdef UWDIAG
         qcu_inv(:idim,k_inv)     = qcu(:idim,k)  ! Diagnostic only
         qlu_inv(:idim,k_inv)     = qlu(:idim,k)
         qiu_inv(:idim,k_inv)     = qiu(:idim,k)
         qc_inv(:idim,k_inv)      = qc(:idim,k)
         xc_inv(:idim,k_inv)      = xc(:idim,k)
#endif
         if (dotransport.eq.1) then
         do m = 1, ncnst
            do i=1,idim
               tr0(i,k,m)   = MAX(mintracer,tr0(i,k,m))
            enddo
            CNV_Tracers(m)%Q(:,:,k_inv) = reshape(tr0(:,k,m), (/IM,JM/))
#ifdef UWDIAG
            trten_inv(:idim,k_inv,m) = trten(:idim,k,m)            
#endif
         enddo
         endif
      end do
      dcm_inv(:idim,k0) = 0.
      ! Re-scale liquid/ice water sub-tendencies to enforce conservation
      where(ABS(qldet_inv+qlsub_inv).gt.1e-12)
        tmp2d = qlten_inv / (qldet_inv+qlsub_inv)
        qldet_inv = tmp2d*qldet_inv
        qlsub_inv = tmp2d*qlsub_inv
      end where
      where(ABS(qidet_inv+qisub_inv).gt.1e-12)
        tmp2d = qiten_inv / (qidet_inv+qisub_inv)
        qidet_inv = tmp2d*qidet_inv
        qisub_inv = tmp2d*qisub_inv
      end where
      

   end subroutine compute_uwshcu_inv


   subroutine compute_uwshcu(idim, k0, dt,ncnst, pifc0_in,zifc0_in,& ! IN
         exnifc0_in, pmid0_in, zmid0_in, exnmid0_in, dp0_in,       &
         u0_in, v0_in, qv0_in, ql0_in, qi0_in, th0_in,             &
         tr0_inout, kpbl_in, tke_in, cush_inout,                   & ! OUT
         umf_out, dcm_out, qvten_out, qlten_out, qiten_out,        &
         sten_out, uten_out, vten_out, qrten_out,                  &
         qsten_out, cufrc_out, fer_out, fdr_out, qldet_out,        &
         qidet_out, qlsub_out, qisub_out, ndrop_out, nice_out,     &
         shfx, evap, cnvtr, tpert_out, qpert_out,                  &
         qtflx_out, slflx_out, uflx_out, vflx_out,                 &
#ifdef UWDIAG
         qcu_out, qlu_out, qiu_out, cbmf_out, qc_out,              & ! DIAG ONLY
         cnt_out, cnb_out, cinh_out, plcl_out, plfc_out, pinv_out, &
         prel_out, pbup_out, wlcl_out, qtsrc_out, thlsrc_out,      &
         thvlsrc_out, tkeavg_out, cldhgt_out, wu_out, qtu_out,     &
         thlu_out, thvu_out, uu_out, vu_out, xc_out, trten_out,    &
#endif
         dotransport)  

    ! ------------------------------------------------------------ !
    !                                                              !  
    !  University of Washington Shallow Convection Scheme          !
    !                                                              !
    !  Described in Park and Bretherton. 2008. J. Climate :        !
    !                                                              !
    ! 'The University of Washington shallow convection and         !
    !  moist turbulent schemes and their impact on climate         !
    !  simulations with the Community Atmosphere Model'            !
    !                                                              !
    !  Coded in CESM by Sungsu Park. Oct.2005.                     ! 
    !                                May.2008.                     !
    !                                                              !
    !  Coded in GEOS by Nathan Arnold. July 2016.                  !
    !                                                              !
    !                                                              !
    !  For general questions, email sungsup@ucar.edu or            ! 
    !                               sungsu@atmos.washington.edu    !
    !                                                              !
    !  For GEOS-specific questions, email nathan.arnold@nasa.gov   !
    !                                                              !
    ! ------------------------------------------------------------ !

      integer, intent(in)  :: idim               ! Number of columns
      integer, intent(in)  :: k0                 ! Number of vertical levels
      integer, intent(in)  :: ncnst              ! Number of tracers
      integer, intent(in)  :: dotransport        ! Transport tracers [1 true]
      real,    intent(in)  :: dt                 ! Timestep [s]

      real, intent(in)    :: pifc0_in(   idim,0:k0 )  ! Environmental pressure at interfaces [Pa]
      real, intent(in)    :: zifc0_in(   idim,0:k0 )  ! Environmental height at interfaces [m]
      real, intent(in)    :: exnifc0_in( idim,0:k0 )  ! Exner function at interfaces
      real, intent(in)    :: pmid0_in(   idim,k0   )  ! Environmental pressure at midpoints [Pa]
      real, intent(in)    :: zmid0_in(   idim,k0   )  ! Environmental height at midpoints [m]
      real, intent(in)    :: exnmid0_in( idim,k0   )  ! Exner function at midpoints
      real, intent(in)    :: dp0_in( idim,k0 )        ! Environmental layer pressure thickness
      real, intent(in)    :: u0_in ( idim,k0 )        ! Environmental zonal wind [m/s]
      real, intent(in)    :: v0_in ( idim,k0 )        ! Environmental meridional wind [m/s]
      real, intent(in)    :: qv0_in( idim,k0 )        ! Environmental specific humidity
      real, intent(in)    :: ql0_in( idim,k0 )        ! Environmental liquid water specific humidity
      real, intent(in)    :: qi0_in( idim,k0 )        ! Environmental ice specific humidity
      real, intent(in)    :: th0_in ( idim,k0 )       ! Environmental potential temperature [K]
      real, intent(in)    :: tke_in( idim,0:k0 )      ! Turbulent kinetic energy at interfaces
      real, intent(in)    :: shfx(idim)               ! Surface sensible heat
      real, intent(in)    :: evap(idim)               ! Surface evaporation
      real, intent(in)    :: cnvtr(idim)              ! Convective tracer
      real, intent(out)   :: tpert_out(idim)          ! Temperature perturbation
      real, intent(out)   :: qpert_out(idim)          ! Humidity perturbation
      real, intent(out)   :: qtflx_out(idim, 0:k0 )   
      real, intent(out)   :: slflx_out(idim, 0:k0 )   
      real, intent(out)   :: uflx_out(idim, 0:k0 )   
      real, intent(out)   :: vflx_out(idim, 0:k0 )   
      integer, intent(in) :: kpbl_in( idim )          ! Boundary layer top layer index

      real, intent(inout) :: cush_inout( idim )       ! Convective scale height [m]
      real, intent(inout) :: tr0_inout(idim,k0,ncnst) !  Environmental tracers [ #, kg/kg ]

      real, intent(out)   :: umf_out(idim,0:k0)       !  Updraft mass flux at the interfaces [ kg/m2/s ]
      real, intent(out)   :: dcm_out(idim,k0)         !  Detrained cloudy air mass
      real, intent(out)   :: qvten_out(idim,k0)       !  Tendency of water vapor specific humidity [ kg/kg/s ]
      real, intent(out)   :: qlten_out(idim,k0)       !  Tendency of liquid water specific humidity [ kg/kg/s ]
      real, intent(out)   :: qiten_out(idim,k0)       !  Tendency of ice specific humidity [ kg/kg/s ]
      real, intent(out)   :: sten_out(idim,k0)        !  Tendency of dry static energy [ J/kg/s ]
      real, intent(out)   :: uten_out(idim,k0)        !  Tendency of zonal wind [ m/s2 ]
      real, intent(out)   :: vten_out(idim,k0)        !  Tendency of meridional wind [ m/s2 ]
      real, intent(out)   :: qrten_out(idim,k0)       !  Tendency of rain water specific humidity [ kg/kg/s ]
      real, intent(out)   :: qsten_out(idim,k0)       !  Tendency of snow specific humidity [ kg/kg/s ]
      real, intent(out)   :: cufrc_out(idim,k0)       !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
      real, intent(out)   :: fer_out(idim,k0)         !  Fractional lateral entrainment rate [ 1/Pa ]
      real, intent(out)   :: fdr_out(idim,k0)         !  Fractional lateral detrainment rate [ 1/Pa ]

      real, intent(out)   :: qldet_out(idim,k0)
      real, intent(out)   :: qidet_out(idim,k0)
      real, intent(out)   :: qlsub_out(idim,k0)
      real, intent(out)   :: qisub_out(idim,k0)
      real, intent(out)   :: ndrop_out(idim,k0)
      real, intent(out)   :: nice_out(idim,k0)

!--------- Diagnostic only ------------
#ifdef UWDIAG
      real, intent(out)   :: trten_out(idim,k0,ncnst) !  Tendency of tracers [ #/s, kg/kg/s ]
      real, intent(out)   :: wu_out(idim,0:k0)        !  Updraft vertical velocity
      real, intent(out)   :: qtu_out(idim,0:k0)       !  Updraft qt [ kg/kg ]
      real, intent(out)   :: thlu_out(idim,0:k0)      !  Updraft thl [ K ]
      real, intent(out)   :: thvu_out(idim,0:k0)      !  Updraft thv [ K ]
      real, intent(out)   :: uu_out(idim,0:k0)        !  Updraft zonal wind [ m/s ] 
      real, intent(out)   :: vu_out(idim,0:k0)        !  Updraft meridional wind [ m/s ]
      real, intent(out)   :: qcu_out(idim,k0)         !  Condensate water specific humidity within cumulus updraft [ kg/kg ]
      real, intent(out)   :: qlu_out(idim,k0)         !  Liquid water specific humidity within cumulus updraft [ kg/kg ]
      real, intent(out)   :: qiu_out(idim,k0)         !  Ice specific humidity within cumulus updraft [ kg/kg ]
      real, intent(out)   :: cbmf_out(idim)           ! Cloud base mass flux [kg/m2/s]
      real, intent(out)   :: qc_out(idim,k0)          !  Tendency of detrained cumulus condensate
      real, intent(out)   :: cnt_out(idim)            ! Cumulus top interface index
      real, intent(out)   :: cnb_out(idim)            ! Cumulus base interface index
      real, intent(out)   :: cinh_out(idim)
      real, intent(out)   :: pinv_out(idim)           !  PBL top pressure [ Pa ]
      real, intent(out)   :: plfc_out(idim)           !  LFC of source air [ Pa ]
      real, intent(out)   :: plcl_out(idim)           !  LCL of source air [ Pa ]
      real, intent(out)   :: prel_out(idim)
      real, intent(out)   :: pbup_out(idim)
      real, intent(out)   :: tkeavg_out(idim)         !  Average tke over the PBL [ m2/s2 ]
      real, intent(out)   :: cldhgt_out(idim)
      real, intent(out)   :: xc_out(idim,k0)
#endif

    !
    ! Internal Output Variables
    !
!    real          qtten_out(idim,k0)             !  Tendency of qt [ kg/kg/s ]
!    real          slten_out(idim,k0)             !  Tendency of sl [ J/kg/s ]
!    real          ufrc_out(idim,0:k0)            !  Updraft fractional area at the interfaces [ fraction ]


      
      !-----------------------------------------------
      ! One-dimensional variables at each grid point
      !-----------------------------------------------

      ! Input variables
 
      real :: pifc0(0:k0)
      real :: zifc0(0:k0)
      real :: pmid0(k0)
      real :: zmid0(k0)
      real :: dp0(k0)
      real :: u0(k0)
      real :: v0(k0)
      real :: tke(1:k0)
      real :: qv0(k0)
      real :: ql0(k0)
      real :: qi0(k0)
      real :: cush
      real :: tr0(k0,ncnst)

      ! Environmental variables derived from input variables

      real :: t0(k0)
      real :: s0(k0)
      real :: qt0(k0)
      real :: thl0(k0)
      real :: thvl0(k0)
      real :: ssqt0(k0)
      real :: ssthl0(k0)
      real :: ssu0(k0)
      real :: ssv0(k0)
      real :: thv0bot(k0)
      real :: thv0top(k0)
      real :: thvl0bot(k0)
      real :: thvl0top(k0)
      real :: exnmid0(k0)
      real :: exnifc0(0:k0)
      real :: sstr0(k0,ncnst)

   ! 2-1. For preventing negative condensate at the provisional time step

    real    qv0_star(k0)                                 !  Environmental water vapor specific humidity [ kg/kg ]
    real    ql0_star(k0)                                 !  Environmental liquid water specific humidity [ kg/kg ]
    real    qi0_star(k0)                                 !  Environmental ice specific humidity [ kg/kg ]
    real    s0_star(k0)                                  !  Environmental dry static energy [ J/kg ]


   ! 3. Variables associated with cumulus convection

    real    umf(0:k0)                                    !  Updraft mass
    real    emf(0:k0)                                    !  Penetrative
    real    dcm(k0)                                      !  Detrained cloud mass
    real    qvten(k0)                                    !  Tendency of
    real    qlten(k0)                                    !  Tendency of
    real    qiten(k0)                                    !  Tendency of
    real    sten(k0)                                     !  Tendency of
    real    uten(k0)                                     !  Tendency of
    real    vten(k0)                                     !  Tendency of
    real    qrten(k0)                                    !  Tendency of
    real    qsten(k0)                                    !  Tendency of
    real    slflx(0:k0)                                  !

    real    qtflx(0:k0)                                  !
    real    uflx(0:k0)                                   !
    real    vflx(0:k0)                                   !
    real    cufrc(k0)                                    !  Shallow
    real    qcu(k0)                                      !  Condensate

    real    qlu(k0)                                      !  Liquid water
    real    qiu(k0)                                      !  Ice specific
    real    dwten(k0)                                    !  Detrained
    real    diten(k0)                                    !  Detrained ice
    real    fer(k0)                                      !  Fractional
    real    fdr(k0)                                      !  Fractional
    real    xco(k0)
    real    uf(k0)                                       !  Zonal wind at
    real    vf(k0)                                       !  Meridional
    real    qc(k0)                                       !  Tendency due
    real    qlten_det(k0)
    real    qiten_det(k0)
!    real    qldet(k0)
!    real    qidet(k0)

    real    qc_l(k0)                                     !  Tendency due
    real    qc_i(k0)                                     !  Tendency due
    real    qc_lm
    real    qc_im
    real    nc_lm
    real    nc_im
    real    totsink
    real    ql_emf_kbup
    real    qi_emf_kbup
    real    nl_emf_kbup
    real    ni_emf_kbup
    real    cnt                                           !  Cumulus top
    real    cnb                                           !  Cumulus base
    real    qtten(k0)                                    !  Tendency of
    real    slten(k0)                                    !  Tendency of
    real    ufrc(0:k0)                                   !  Updraft
    real    trten(k0,ncnst)                              !  Tendency of
    real    trflx(0:k0,ncnst)                            !  Flux of
    real    trflx_d(0:k0)                                !  Adjustive
    real    trflx_u(0:k0)                                !  Adjustive
    real    trmin
    real    pdelx, dum 

! Variables for temperature/moisture excess in source parcel
    real    zrho, delzg,buoyflx,wstar !pahfs,pqhfl,zkhvfl,pgeoh,zws,ztexec,zqexec

    !----- Variables used for the calculation of condensation sink associated with compensating subsidence
    !      In the current code, this 'sink' tendency is simply set to be zero.

    real    uemf(0:k0)                                   !  Net updraft mass flux at the interface ( emf + umf ) [ kg/m2/s ]
    real    comsub(k0)                                   !  Compensating subsidence
                                                              ! at the layer mid-point ( unit of mass flux, umf ) [ kg/m2/s ]
    real    qlten_sink(k0)                               !  Liquid condensate tendency
                                                              ! by compensating subsidence/upwelling [ kg/kg/s ]
    real    qiten_sink(k0)                               !  Ice    condensate tendency
                                                              ! by compensating subsidence/upwelling [ kg/kg/s ]
    real    thlten_sub, qtten_sub                         !  Tendency of conservative scalars
                                                              ! by compensating subsidence/upwelling
    real    qlten_sub, qiten_sub                          !  Tendency of ql0, qi0
                                                              ! by compensating subsidence/upwelling
    real    nlten_sub, niten_sub                          !  Tendency of nl0, ni0
                                                              ! by compensating subsidence/upwelling
    real    thl_prog, qt_prog                             !  Prognosed 'thl, qt'
                                                              ! by compensating subsidence/upwelling


    !----- Variables describing cumulus updraft

    real    wu(0:k0)                                     !  Updraft vertical velocity at the interface [ m/s ]
    real    thlu(0:k0)                                   !  Updraft liquid potential temperature at the interface [ K ]
    real    qtu(0:k0)                                    !  Updraft total specific humidity at the interface [ kg/kg ]
    real    uu(0:k0)                                     !  Updraft zonal wind at the interface [ m/s ]
    real    vu(0:k0)                                     !  Updraft meridional wind at the interface [ m/s ]
    real    thvu(0:k0)                                   !  Updraft virtual potential temperature at the interface [ m/s ]
    real    rei(k0)                                      !  Updraft fractional mixing rate with the environment [ 1/Pa ]
    real    tru(0:k0,ncnst)                              !  Updraft tracers [ #, kg/kg ]

    !----- Variables describing conservative scalars of entraining downdrafts  at the 
    !      entraining interfaces, i.e., 'kbup <= k < kpen-1'. At the other interfaces,
    !      belows are simply set to equal to those of updraft for simplicity - but it
    !      does not influence numerical calculation.

    real    thlu_emf(0:k0)                               !  Penetrative downdraft liquid potential temperature
                                                              ! at entraining interfaces [ K ]
    real    qtu_emf(0:k0)                                !  Penetrative downdraft total water
                                                              ! at entraining interfaces [ kg/kg ]
    real    uu_emf(0:k0)                                 !  Penetrative downdraft zonal wind
                                                              ! at entraining interfaces [ m/s ]
    real    vu_emf(0:k0)                                 !  Penetrative downdraft meridional wind
                                                              ! at entraining interfaces [ m/s ]
    real    tru_emf(0:k0,ncnst)                          !  Penetrative Downdraft tracers
                                                              ! at entraining interfaces [ #, kg/kg ] 


      ! Other internal variables

      integer    kk, k, i, kp1, km1, mm, m
      integer    iter_scaleh, iter_xc
      integer    id_check, status

      integer    klcl        ! Layer containing LCL of source air
      integer    kinv        ! Inversion layer with PBL top interface as lower interface
      integer    krel        ! Release layer where buoyancy sorting first occurs
      integer    klfc        ! LFC layer of cumulus source air
      integer    kbup        ! Top layer in which buoyancy is positive at top interface
      integer    kpen        ! Highest layer with positive updraft velocity

      logical    id_exit
      logical    forcedCu

      real       cin, cinlcl
      real       thlsrc, qtsrc, usrc, vsrc, thvlsrc
      real       uplus, vplus
      real       trsrc(ncnst), tre(ncnst)
      real       plcl, plfc, prel, wrel
      real       ee2, ud2, wtw, wtwb
      real       xc
      real       cldhgt, scaleh, tscaleh, cridis
      real       sigmaw, tkeavg, qtavg, thvlavg, uavg, vavg, dpsum, dpi, thvlmin
      real       thlxsat, qtxsat, thvxsat, x_cu, x_en, thv_x0, thv_x1
      real       dpe, exne, thvebot, thle, qte, ue, ve, thlue, qtue, wue
      real       mu, mumin0, mumin2, mulcl, mulclstar
      real       cbmf, wcrit, winv, wlcl, ufrcinv, ufrclcl
      real       exql, exqi, ppen

      real       thj, qvj, qlj, qij, thvj, tj, thv0j, rhomid0j, rhoifc0j, qse
      real       thl0top, thl0bot, qt0bot, qt0top, thvubot, thvutop
      real       thl0lcl, qt0lcl, thv0lcl, thv0rel, rho0inv, autodet
      real       thlu_top, qtu_top, qlu_top, qiu_top, qlu_mid, qiu_mid, exntop
      real       aquad, bquad, cquad, xc1, xc2, excessu, excess0, xsat, xs1, xs2
      real       bogbot, bogtop, delbog, drage, expfac
      real       rcwp, rlwp, riwp, qcubelow, qlubelow, qiubelow
      real       qs
      real         qsat_arg, qsat_pe 
      real       pe
      real       xsrc, xmean, xtop, xbot, xflx(0:k0)

     
    !----- Some diagnostic internal output variables

#ifdef UWDIAG
    real  trflx_out(idim,0:k0,ncnst)           !  Updraft/pen.entrainment tracer flux [ #/m2/s, kg/kg/m2/s ] 
    real  ufrcinvbase_out(idim)                !  Cumulus updraft fraction at the PBL top [ fraction ]
    real  ufrclcl_out(idim)                    !  Cumulus updraft fraction at the LCL
                                               ! ( or PBL top when LCL is below PBL top ) [ fraction ]
    real  winvbase_out(idim)                   !  Cumulus updraft velocity at the PBL top [ m/s ]
    real  wlcl_out(idim)                       !  Cumulus updraft velocity at the LCL
                                               ! ( or PBL top when LCL is below PBL top ) [ m/s ]
!    real  pbup_out(idim)                      !  Highest interface level of positive buoyancy [ Pa ]
    real  ppen_out(idim)                       !  Highest interface evel where Cu w = 0 [ Pa ]
    real  qtsrc_out(idim)                      !  Source air qt [ kg/kg ]
    real  thlsrc_out(idim)                     !  Source air thl [ K ]
    real  thvlsrc_out(idim)                    !  Source air thvl [ K ]
    real  emfkbup_out(idim)                    !  Penetrative downward mass flux at 'kbup' interface [ kg/m2/s ]
    real  cinlclh_out(idim)                    !  Convective INhibition upto LCL (CIN) [ J/kg = m2/s2 ]
    real  cbmflimit_out(idim)                  !  Cloud base mass flux limiter [ kg/m2/s ]
    real  zinv_out(idim)                       !  PBL top height [ m ]
    real  rcwp_out(idim)                       !  Layer mean Cumulus LWP+IWP [ kg/m2 ] 
    real  rlwp_out(idim)                       !  Layer mean Cumulus LWP [ kg/m2 ] 
    real  riwp_out(idim)                       !  Layer mean Cumulus IWP [ kg/m2 ] 

    real  qtu_emf_out(idim,0:k0)               !  Penetratively entrained qt [ kg/kg ]   
    real  thlu_emf_out(idim,0:k0)              !  Penetratively entrained thl [ K ]
    real  uu_emf_out(idim,0:k0)                !  Penetratively entrained u [ m/s ]
    real  vu_emf_out(idim,0:k0)                !  Penetratively entrained v [ m/s ]
    real  uemf_out(idim,0:k0)                  !  Net upward mass flux
                                               ! including penetrative entrainment (umf+emf) [ kg/m2/s ]
    real  dwten_out(idim,k0)
    real  diten_out(idim,k0)
    real  tru_out(idim,0:k0,ncnst)             !  Updraft tracers [ #, kg/kg ]   
    real  tru_emf_out(idim,0:k0,ncnst)         !  Penetratively entrained tracers [ #, kg/kg ]
    real  wu_s(0:k0)                           !  Same as above but for implicit CIN
    real  qtu_s(0:k0)
    real  thlu_s(0:k0)
    real  thvu_s(0:k0)
    real  uu_s(0:k0)
    real  vu_s(0:k0)
    real  qtu_emf_s(0:k0) 
    real  thlu_emf_s(0:k0)  
    real  uu_emf_s(0:k0)   
    real  vu_emf_s(0:k0)
    real  uemf_s(0:k0)   
    real  tru_s(0:k0,ncnst)
    real  tru_emf_s(0:k0,ncnst)   

    real  dwten_s(k0)
    real  diten_s(k0)

    real  excessu_arr_out(idim,k0)
    real  excessu_arr(k0) 
    real  excessu_arr_s(k0)
    real  excess0_arr_out(idim,k0)
    real  excess0_arr(k0)
    real  excess0_arr_s(k0)
    real  xc_arr_out(idim,k0)
    real  xc_arr(k0)
    real  xc_arr_s(k0)
    real  aquad_arr_out(idim,k0)
    real  aquad_arr(k0)
    real  aquad_arr_s(k0)
    real  bquad_arr_out(idim,k0)
    real  bquad_arr(k0)
    real  bquad_arr_s(k0)
    real  cquad_arr_out(idim,k0) 
    real  cquad_arr(k0)
    real  cquad_arr_s(k0)
    real  bogbot_arr_out(idim,k0)
    real  bogbot_arr(k0)
    real  bogbot_arr_s(k0)
    real  bogtop_arr_out(idim,k0)
    real  bogtop_arr(k0)
    real  bogtop_arr_s(k0)
#endif

    real       exit_ufrc(idim)
    real       exit_wtw(idim)
    real       exit_drycore(idim)
    real       exit_wu(idim)
    real       exit_cufilter(idim)
    real       exit_rei(idim)
    real       exit_kinv1(idim)
    real       exit_klfck0(idim)
    real       exit_klclk0(idim)
    real       exit_uwcu(idim)
    real       exit_conden(idim)

    real       limit_cinlcl(idim)
    real       limit_cin(idim)
    real       ind_delcin(idim)
    real       limit_rei(idim)
    real       limit_shcu(idim)
    real       limit_negcon(idim)
    real       limit_ufrc(idim)
    real       limit_ppen(idim)
    real       limit_emf(idim)
    real       limit_cbmf(idim)

    real :: ufrcinvbase_s, ufrclcl_s, winvbase_s, wlcl_s, plcl_s, pinv_s, prel_s, plfc_s, &
              qtsrc_s, thlsrc_s, thvlsrc_s, emfkbup_s, cinlcl_s, pbup_s, ppen_s, cbmflimit_s, &
              tkeavg_s, zinv_s, rcwp_s, rlwp_s, riwp_s 
    real :: ufrcinvbase, winvbase, emfkbup, cbmflimit

    !----- Variables for implicit CIN computation

    real, dimension(k0)   :: dcm_s, qv0_s  , ql0_s   , qi0_s   , s0_s    , u0_s    ,          & 
                             v0_s   , t0_s    , & !qt0_s   , thl0_s  ,thvl0_s ,           &
                             qvten_s , qlten_s, qiten_s , qrten_s , qsten_s , sten_s  ,& 
                             uten_s, vten_s, cufrc_s, qcu_s, qlu_s, qiu_s, & 
                             fer_s, fdr_s, xc_s, qc_s, qtten_s, slten_s, qldet_s,        &
                             qidet_s, qlsub_s, qisub_s
    real, dimension(0:k0) :: umf_s  , slflx_s , qtflx_s , ufrc_s  , uflx_s , vflx_s
    real                  :: cush_s , cin_s   , cbmf_s, cnt_s, cnb_s
    real                  :: cin_i, cin_f, del_CIN, ke, alpha !, thlj
    real                  :: cinlcl_i, cinlcl_f, del_cinlcl
    integer                 :: iter
    
    real, dimension(k0,ncnst)   :: tr0_s, trten_s
    real, dimension(0:k0,ncnst) :: trflx_s

    !----- Variables for temporary storage

      real, dimension(k0) :: qv0_o, ql0_o, qi0_o, t0_o, s0_o, u0_o, v0_o
      real, dimension(k0) :: qt0_o    , thl0_o   , thvl0_o   ,                         &
                              thv0bot_o, thv0top_o, thvl0bot_o, thvl0top_o,             &
                              ssthl0_o , ssqt0_o  , ssu0_o    , ssv0_o
    real                  :: tkeavg_o , thvlmin_o, qtsrc_o, thvlsrc_o, thlsrc_o ,    &
                               usrc_o, vsrc_o, plcl_o, plfc_o, thv0lcl_o
                                        
    integer               :: kinv_o   , klcl_o   , klfc_o  

    real                  :: lts

    real, dimension(k0,ncnst)   :: tr0_o
    real, dimension(k0,ncnst)   :: trten_o, sstr0_o  
    real, dimension(0:k0,ncnst) :: trflx_o
    real, dimension(ncnst)       :: trsrc_o
    integer                          :: ixnumliq, ixnumice, ixcldliq, ixcldice


    ! ------------------ !
    !                    !
    ! Define Parameters  !
    !                    !
    ! ------------------ !

    ! ------------------------ !
    ! Source air perturbations !
    ! ------------------------ !

    real :: qtsrc_fac 
    real :: thlsrc_fac
    real :: qtsrchgt

    ! ------------------------ !
    ! Iterative xc calculation !
    ! ------------------------ !

    integer :: niter_xc   

    ! ----------------------------------------------------------- !
    ! Choice of 'CIN = cin' (.true.) or 'CIN = cinlcl' (.false.). !
    ! ----------------------------------------------------------- !

    logical :: use_CINcin

    ! --------------------------------------------------------------- !
    ! Choice of 'explicit' ( 1 ) or 'implicit' ( 2 )  CIN.            !
    !                                                                 !
    ! When choose 'CIN = cinlcl' above,  it is recommended not to use ! 
    ! implicit CIN, i.e., do 'NOT' choose simultaneously :            !
    !            [ 'use_CINcin=.false. & 'iter_cin=2' ]               !
    ! since 'cinlcl' will be always set to zero whenever LCL is below !
    ! the PBL top interface in the current code. So, averaging cinlcl !
    ! of two iter_cin steps is likely not so good. Except that,   all !
    ! the other combinations of  'use_CINcin'  & 'iter_cin' are OK.   !
    ! --------------------------------------------------------------- !

    integer :: iter_cin 

    ! ---------------------------------------------------------------- !
    ! Choice of 'self-detrainment' by negative buoyancy in calculating !
    ! cumulus updraft mass flux at the top interface in each layer.    !
    ! ---------------------------------------------------------------- !

    logical :: use_self_detrain
    
    ! --------------------------------------------------------- !
    ! Cumulus momentum flux : turn-on (.true.) or off (.false.) !
    ! --------------------------------------------------------- !

    logical :: use_momenflx

    ! ----------------------------------------------------------------------------------------- !
    ! Penetrative Entrainment : Cumulative ( .true. , original ) or Non-Cumulative ( .false. )  !
    ! This option ( .false. ) is designed to reduce the sensitivity to the vertical resolution. !
    ! ----------------------------------------------------------------------------------------- !

    logical :: use_cumpenent
    logical :: scverbose
    logical :: windsrcavg
    real    :: rpen     ! penetrative entrainment efficiency


    integer :: cridist_opt

    ! ----------------------- !
    ! For lateral entrainment !
    ! ----------------------- !

    real :: rle          !  For critical stopping distance for lateral entrainment [no unit]
    real :: rkm          !  Determine the amount of air that is involved in buoyancy-sorting [no unit]
    real :: mixscale     !  Specify vertical structure of mixing rate
    real :: detrhgt      !  Mixing rate increases above this height to speed detrainment
    real :: rkfre        !  Vertical velocity variance as fraction of  tke. 
    real :: rmaxfrac     !  Maximum allowable 'core' updraft fraction
    real :: mumin1       !  Normalized CIN ('mu') corresponding to 'rmaxfrac' at the PBL top
                         !  obtaind by inverting 'rmaxfrac = 0.5*erfc(mumin1)'.
                         !  [rmaxfrac:mumin1]=[ 0.05:1.163, 0.075:1.018, 0.1:0.906, 0.15:0.733, 0.2:0.595, 0.25:0.477]
    real :: rbuoy        !  For nonhydrostatic pressure effects on updraft [no unit]
    real :: rdrag        !  Drag coefficient [no unit]
    real :: epsvarw      !  Variance of w at PBL top by meso-scale component [m2/s2]          
    real :: PGFc         !  This is used for calculating vertical variations cumulus  
                         !  'u' & 'v' by horizontal PGF during upward motion [no unit]
    real :: frc_rasn

!!! TEMPORARY:  should be ncnst array of minimum values for all constituents
    real, parameter,dimension(4) :: qmin = [0.,0.,0.,0.]


    ! ---------------------------------------- !
    ! Bulk microphysics controlling parameters !
    ! --------------------------------------------------------------------------- ! 
    ! criqc    : Maximum condensate that can be hold by cumulus updraft [kg/kg]   !
    ! frc_rasn : Fraction of precipitable condensate in the expelled cloud water  !
    !            from cumulus updraft. The remaining fraction ('1-frc_rasn')  is  !
    !            'suspended condensate'.                                          !
    !                0 : all expelled condensate is 'suspended condensate'        ! 
    !                1 : all expelled condensate is 'precipitable condensate'     !
    ! kevp     : Evaporative efficiency                                           !
    ! noevap_krelkpen : No evaporation from 'krel' to 'kpen' layers               ! 
    ! --------------------------------------------------------------------------- !    

    real :: criqc
    real :: rdrop

   ! ------------------------ !
   ! Assign parameter values  !
   ! ------------------------ !

    niter_xc         = shlwparams%niter_xc
    use_CINcin       = shlwparams%use_CINcin
    iter_cin         = shlwparams%iter_cin
    use_self_detrain = shlwparams%use_self_detrain
    use_momenflx     = shlwparams%use_momenflx
    use_cumpenent    = shlwparams%use_cumpenent
    scverbose        = shlwparams%scverbose
    windsrcavg       = shlwparams%windsrcavg
    rpen             = shlwparams%rpen
    cridist_opt      = shlwparams%cridist_opt
    rle       = shlwparams%rle      !  For critical stopping distance for lateral entrainment [no unit]
    rkm       = shlwparams%rkm      !  Determine the amount of air that is involved in buoyancy-sorting [no unit]
    mixscale  = shlwparams%mixscale !  Specifies vertical structure of mixing rate
    detrhgt   = shlwparams%detrhgt  !  Specifies vertical structure of mixing rate
    rkfre     = shlwparams%rkfre    !  Vertical velocity variance as fraction of  tke. 
    rmaxfrac  = shlwparams%rmaxfrac !  Maximum allowable 'core' updraft fraction
    mumin1    = shlwparams%mumin1
    rbuoy     = shlwparams%rbuoy    !  For nonhydrostatic pressure effects on updraft [no unit]
    rdrag     = shlwparams%rdrag    !  Drag coefficient [no unit]
    epsvarw   = shlwparams%epsvarw  !  Variance of w at PBL top by meso-scale component [m2/s2]          
    PGFc      = shlwparams%PGFc     !  This is used for calculating vertical variations cumulus  
                                    !  'u' & 'v' by horizontal PGF during upward motion [no unit]
    criqc     = shlwparams%criqc
    frc_rasn  = shlwparams%frc_rasn
    rdrop     = shlwparams%rdrop
    thlsrc_fac = shlwparams%thlsrc_fac
    qtsrc_fac  = shlwparams%qtsrc_fac
    qtsrchgt   = shlwparams%qtsrchgt

    !------------------------!
    !                        !
    ! Start Main Calculation !
    !                        !
    !------------------------!

!!! TEMPORARY test settings
    ixcldice = 1
    ixcldliq = 2
    ixnumliq = 3
    ixnumice = 4


    ! ------------------------------------------------------- !
    ! Initialize output variables defined for all grid points !
    ! ------------------------------------------------------- !

    umf_out(:idim,0:k0)          = 0.0
    dcm_out(:idim,:k0)           = 0.0
    cufrc_out(:idim,:k0)         = 0.0
    fer_out(:idim,:k0)           = MAPL_UNDEF
    fdr_out(:idim,:k0)           = MAPL_UNDEF
    qldet_out(:idim,:k0)         = 0.0
    qidet_out(:idim,:k0)         = 0.0
    qlsub_out(:idim,:k0)         = 0.0
    qisub_out(:idim,:k0)         = 0.0
    ndrop_out(:idim,:k0)         = 0.0
    nice_out(:idim,:k0)          = 0.0
    qtflx_out(:idim,0:k0)        = 0.0
    slflx_out(:idim,0:k0)        = 0.0
    uflx_out(:idim,0:k0)         = 0.0
    vflx_out(:idim,0:k0)         = 0.0
    tpert_out(:idim)             = 0.0
    qpert_out(:idim)             = 0.0

#ifdef UWDIAG
    cbmf_out(:idim)              = 0.0
    cinh_out(:idim)              = MAPL_UNDEF
    cinlclh_out(:idim)           = MAPL_UNDEF
    cldhgt_out(:idim)            = 0.0
    qcu_out(:idim,:k0)           = 0.0
    qlu_out(:idim,:k0)           = 0.0
    qiu_out(:idim,:k0)           = 0.0
    qc_out(:idim,:k0)            = 0.0
    cnt_out(:idim)               = real(k0)
    cnb_out(:idim)               = 0.0
    xc_out(:idim,:k0)            = 0.0
!    ufrc_out(:idim,0:k0)         = 0.0
!    uflx_out(:idim,0:k0)         = 0.0
!    vflx_out(:idim,0:k0)         = 0.0

    ufrcinvbase_out(:idim)       = 0.0
    ufrclcl_out(:idim)           = 0.0
    winvbase_out(:idim)          = 0.0
    wlcl_out(:idim)              = 0.0
    plcl_out(:idim)              = 0.0
    pinv_out(:idim)              = 0.0
    plfc_out(:idim)              = 0.0
    prel_out(:idim)              = 0.0
    pbup_out(:idim)              = 0.0
    ppen_out(:idim)              = 0.0
    qtsrc_out(:idim)             = 0.0
    thlsrc_out(:idim)            = 0.0
    thvlsrc_out(:idim)           = 0.0
    emfkbup_out(:idim)           = 0.0
    cbmflimit_out(:idim)         = 0.0
    tkeavg_out(:idim)            = 0.0
    zinv_out(:idim)              = 0.0
    rcwp_out(:idim)              = 0.0
    rlwp_out(:idim)              = 0.0
    riwp_out(:idim)              = 0.0
  
    wu_out(:idim,0:k0)          = MAPL_UNDEF
    qtu_out(:idim,0:k0)         = MAPL_UNDEF
    thlu_out(:idim,0:k0)        = MAPL_UNDEF
    thvu_out(:idim,0:k0)        = MAPL_UNDEF
    uu_out(:idim,0:k0)          = MAPL_UNDEF
    vu_out(:idim,0:k0)          = MAPL_UNDEF
    qtu_emf_out(:idim,0:k0)     = 0.0
    thlu_emf_out(:idim,0:k0)    = 0.0
    uu_emf_out(:idim,0:k0)      = 0.0
    vu_emf_out(:idim,0:k0)      = 0.0
    uemf_out(:idim,0:k0)        = 0.0

    dwten_out(:idim,:k0)        = 0.0
    diten_out(:idim,:k0)        = 0.0

    trten_out(:idim,:k0,:ncnst)    = 0.0
    trflx_out(:idim,0:k0,:ncnst)   = 0.0
    tru_out(:idim,0:k0,:ncnst)     = 0.0
    tru_emf_out(:idim,0:k0,:ncnst) = 0.0

    excessu_arr_out(:idim,:k0)   = 0.0
    excess0_arr_out(:idim,:k0)   = 0.0
    xc_arr_out(:idim,:k0)        = 0.0
    aquad_arr_out(:idim,:k0)     = 0.0
    bquad_arr_out(:idim,:k0)     = 0.0
    cquad_arr_out(:idim,:k0)     = 0.0
    bogbot_arr_out(:idim,:k0)    = 0.0
    bogtop_arr_out(:idim,:k0)    = 0.0
#endif

    exit_UWCu(:idim)             = 0.0 
    exit_conden(:idim)           = 0.0 
    exit_klclk0(:idim)           = 0.0 
    exit_klfck0(:idim)           = 0.0 
    exit_ufrc(:idim)             = 0.0 
    exit_wtw(:idim)              = 0.0 
    exit_drycore(:idim)          = 0.0 
    exit_wu(:idim)               = 0.0 
    exit_cufilter(:idim)         = 0.0 
    exit_kinv1(:idim)            = 0.0 
    exit_rei(:idim)              = 0.0 

    limit_shcu(:idim)            = 0.0 
    limit_negcon(:idim)          = 0.0 
    limit_ufrc(:idim)            = 0.0
    limit_ppen(:idim)            = 0.0
    limit_emf(:idim)             = 0.0
    limit_cinlcl(:idim)          = 0.0
    limit_cin(:idim)             = 0.0
    limit_cbmf(:idim)            = 0.0
    limit_rei(:idim)             = 0.0

    ind_delcin(:idim)            = 0.0


      !========================
      !   Start column loop
      !========================

      do i = 1, idim

         id_exit = .false.

         pifc0(0:k0)     = pifc0_in(i,0:k0)
         zifc0(0:k0)     = zifc0_in(i,0:k0)
         pmid0(:k0)      = pmid0_in(i,:k0)
         zmid0(:k0)      = zmid0_in(i,:k0)
         dp0(:k0)        = dp0_in(i,:k0)
         u0(:k0)         = u0_in(i,:k0)
         v0(:k0)         = v0_in(i,:k0)
         qv0(:k0)        = qv0_in(i,:k0)
         ql0(:k0)        = ql0_in(i,:k0)
         qi0(:k0)        = qi0_in(i,:k0)
         tke(1:k0)       = tke_in(i,1:k0)
!         pblh            = pblh_in(i)
         cush            = cush_inout(i)

         if (dotransport.eq.1) then
         do m = 1,ncnst   ! loop over tracers
            tr0(:k0,m) = tr0_inout(i,:k0,m)
         end do
         endif

         !------------------------------------------------------!
         ! Compute basic thermodynamic variables directly from  !
         ! input variables for each column                      !
         !------------------------------------------------------!

         ! Compute internal environmental variables

         exnmid0(:k0) = exnmid0_in(i,:k0)
         exnifc0(:k0) = exnifc0_in(i,:k0)
         t0(:k0)      = th0_in(i,:k0) * exnmid0(:k0)
         s0(:k0)      = g*zmid0(:k0) + cp*t0(:k0)
         qt0(:k0)     = qv0(:k0) + ql0(:k0) + qi0(:k0)
         thl0(:k0)    = ( t0(:k0) - xlv*ql0(:k0)/cp - xls*qi0(:k0)/cp ) / exnmid0(:k0)
         thvl0(:k0)   = ( 1. + zvir*qt0(:k0) )*thl0(:k0)


         ! Compute slopes of environmental variables in each layer

         ssthl0    = slope( k0, thl0, pmid0 )
         ssqt0     = slope( k0, qt0 , pmid0 )
         ssu0      = slope( k0, u0  , pmid0 )
         ssv0      = slope( k0, v0  , pmid0 )
         if (dotransport.eq.1) then
         do m = 1,ncnst    ! loop over tracers
            sstr0(:k0,m) = slope( k0, tr0(:k0,m), pmid0 )
         end do
         end if

         ! Compute thv0 and thvl0 at top/bottom interfaces in each layer

         do k = 1,k0

            thl0bot = thl0(k) + ssthl0(k)*(pifc0(k-1) - pmid0(k))
            qt0bot  = qt0(k) + ssqt0(k)*(pifc0(k-1) - pmid0(k))
            call conden( pifc0(k-1),thl0bot,qt0bot,thj,qvj,qlj,qij,qse,id_check )
            if ( id_check .eq. 1 ) then
               exit_conden(i) = 1.0
               id_exit = .true.
               if (scverbose) then
                  call write_parallel('------- UW ShCu: Exit, conden')
               end if
               go to 333
            end if
            thv0bot(k)  = thj*(1. + zvir*qvj - qlj - qij)
            thvl0bot(k) = thl0bot*(1. + zvir*qt0bot)

            thl0top = thl0(k) + ssthl0(k)*(pifc0(k) - pmid0(k))
            qt0top  = qt0(k) + ssqt0(k)*(pifc0(k) - pmid0(k))
            if (k.lt.k0) then
              call conden( pifc0(k),thl0top,qt0top,thj,qvj,qlj,qij,qse,id_check )
              if ( id_check .eq. 1 ) then
                 exit_conden(i) = 1.0
                 id_exit = .true.
                 if (scverbose) then
                   call write_parallel('------- UW ShCu: Exit, conden')
                 end if
                 go to 333
              end if
              thv0top(k)  = thj*(1. + zvir*qvj - qlj - qij)
              thvl0top(k) = thl0top*(1. + zvir*qt0top)
            else
              thv0top(k)  = thv0bot(k)
              thvl0top(k) = thvl0bot(k)
            end if

         end do

!         print *,'thvl0bot=',thvl0bot
!         print *,'thvl0top=',thvl0top

         ! ------------------------------------------------------------ !
         ! Save input and related environmental thermodynamic variables !
         ! for use at "iter_cin=2" when "del_CIN >= 0"                  !
         ! ------------------------------------------------------------ !

         qv0_o(:k0)          = qv0(:k0)
         ql0_o(:k0)          = ql0(:k0)
         qi0_o(:k0)          = qi0(:k0)
         t0_o(:k0)           = t0(:k0)
         s0_o(:k0)           = s0(:k0)
         u0_o(:k0)           = u0(:k0)
         v0_o(:k0)           = v0(:k0)
         qt0_o(:k0)          = qt0(:k0)
         thl0_o(:k0)         = thl0(:k0)
         thvl0_o(:k0)        = thvl0(:k0)
         ssthl0_o(:k0)       = ssthl0(:k0)
         ssqt0_o(:k0)        = ssqt0(:k0)
         thv0bot_o(:k0)      = thv0bot(:k0)
         thv0top_o(:k0)      = thv0top(:k0)
         thvl0bot_o(:k0)     = thvl0bot(:k0)
         thvl0top_o(:k0)     = thvl0top(:k0)
         ssu0_o(:k0)         = ssu0(:k0) 
         ssv0_o(:k0)         = ssv0(:k0) 
         if (dotransport.eq.1) then
         do m = 1, ncnst
           tr0_o(:k0,m)     = tr0(:k0,m)
           sstr0_o(:k0,m)   = sstr0(:k0,m)
         enddo 
         end if
 
         ! ---------------------------------------------- !
         ! Initialize output variables at each grid point !
         ! ---------------------------------------------- !

      umf(0:k0)          = 0.0
      dcm(:k0)           = 0.0
      emf(0:k0)          = 0.0
      slflx(0:k0)        = 0.0
      qtflx(0:k0)        = 0.0
      uflx(0:k0)         = 0.0
      vflx(0:k0)         = 0.0
      qvten(:k0)         = 0.0
      qlten(:k0)         = 0.0
      qiten(:k0)         = 0.0
      sten(:k0)          = 0.0
      uten(:k0)          = 0.0
      vten(:k0)          = 0.0
      qrten(:k0)         = 0.0
      qsten(:k0)         = 0.0
      dwten(:k0)         = 0.0
      diten(:k0)         = 0.0
      cufrc(:k0)         = 0.0
      qcu(:k0)           = 0.0
      qlu(:k0)           = 0.0
      qiu(:k0)           = 0.0
      fer(:k0)           = 0.0
      fdr(:k0)           = 0.0
      xco(:k0)           = 0.0
      cin                = 0.0
      cinlcl             = 0.0
      cbmf               = 0.0
      qc(:k0)            = 0.0
!      qldet(:k0)         = 0.0
!      qidet(:k0)         = 0.0
      qc_l(:k0)          = 0.0
      qc_i(:k0)          = 0.0
      cnt                = real(k0)
      cnb                = 0.0
      qtten(:k0)         = 0.0
      slten(:k0)         = 0.0   
      ufrc(0:k0)         = 0.0  

      thlu(0:k0)         = MAPL_UNDEF
      qtu(0:k0)          = MAPL_UNDEF
      uu(0:k0)           = MAPL_UNDEF
      vu(0:k0)           = MAPL_UNDEF
      wu(0:k0)           = MAPL_UNDEF
      thvu(0:k0)         = MAPL_UNDEF
      thlu_emf(0:k0)     = MAPL_UNDEF
      qtu_emf(0:k0)      = MAPL_UNDEF
      uu_emf(0:k0)       = MAPL_UNDEF
      vu_emf(0:k0)       = MAPL_UNDEF
      
      ufrcinvbase         = 0.0
      ufrclcl             = 0.0
      winvbase            = 0.0
      wlcl                = 0.0
      emfkbup             = 0.0 
      cbmflimit           = 0.0
#ifdef UWDIAG
      excessu_arr(:k0)   = 0.0
      excess0_arr(:k0)   = 0.0
      xc_arr(:k0)        = 0.0
      aquad_arr(:k0)     = 0.0
      bquad_arr(:k0)     = 0.0
      cquad_arr(:k0)     = 0.0
      bogbot_arr(:k0)    = 0.0
      bogtop_arr(:k0)    = 0.0
#endif

      uemf(0:k0)         = 0.0
      comsub(:k0)        = 0.0
      qlten_sink(:k0)    = 0.0
      qiten_sink(:k0)    = 0.0
      qlten_det(:k0)     = 0.0
      qiten_det(:k0)     = 0.0 
!      nlten_sink(:k0)    = 0.0
!      niten_sink(:k0)    = 0.0 

      if (dotransport.eq.1) then
      do m = 1, ncnst
         trflx(0:k0,m)   = 0.0
         trten(:k0,m)    = 0.0
         tru(0:k0,m)     = 0.0
         tru_emf(0:k0,m) = 0.0
      enddo
      endif


        !-----------------------------------------------! 
        ! Below 'iter' loop is for implicit CIN closure !
        !-----------------------------------------------! 

        ! ----------------------------------------------------------------------------- ! 
        ! It is important to note that this iterative cin loop is located at the outer  !
        ! shell of the code. Thus, source air properties can also be changed during the !
        ! iterative cin calculation, because cumulus convection induces non-zero fluxes !
        ! even at interfaces below PBL top height through 'fluxbelowinv' subroutine.    !
        ! ----------------------------------------------------------------------------- !

        do iter = 1, iter_cin

           ! ---------------------------------------------------------------------- ! 
           ! Cumulus scale height                                                   ! 
           ! In contrast to the premitive code, cumulus scale height is iteratively !
           ! calculated at each time step, and at each iterative cin step.          !
           ! It is not clear whether I should locate below two lines within or  out !
           ! of the iterative cin loop.                                             !
           ! ---------------------------------------------------------------------- !

           tscaleh = cush                        
           cush    = -1.
           tkeavg   = 0.
           qtavg   = 0.
           uavg    = 0.
           vavg    = 0.

           ! ----------------------------------------------------------------------- !
           ! Find PBL top height interface index, 'kinv-1' where 'kinv' is the layer !
           ! index with PBLH in it. When PBLH is exactly at interface, 'kinv' is the !
           ! layer index having PBLH as a lower interface.                           !
           ! In the previous code, I set the lower limit of 'kinv' by 2  in order to !
           ! be consistent with the other parts of the code. However in the modified !
           ! code, I allowed 'kinv' to be 1 & if 'kinv = 1', I just exit the program !
           ! without performing cumulus convection. This new approach seems to be    !
           ! more reasonable: if PBL height is within 'kinv=1' layer, surface is STL !
           ! interface (bflxs <= 0) and interface just above the surface should be   !
           ! either non-turbulent (Ri>0.19) or stably turbulent (0<=Ri<0.19 but this !
           ! interface is identified as a base external interface of upperlying CL.  !
           ! Thus, when 'kinv=1', PBL scheme guarantees 'bflxs <= 0'.  For this case !
           ! it is reasonable to assume that cumulus convection does not happen.     !
           ! When these is SBCL, PBL height from the PBL scheme is likely to be very !
           ! close at 'kinv-1' interface, but not exactly, since 'zi' information is !
           ! changed between two model time steps. In order to ensure correct identi !
           ! fication of 'kinv' for general case including SBCL, I imposed an offset !
           ! of 5 [m] in the below 'kinv' finding block.                             !
           ! ----------------------------------------------------------------------- !

           ! invert kpbl index
           kinv = k0 - kpbl_in(i) + 1

15         continue    

           if( kinv .le. 1 ) then        
              exit_kinv1(i) = 1.
              id_exit = .true.
              if (scverbose) then
                call write_parallel('------- UW ShCu: Exit, kinv<=1')
              end if
              go to 333
           endif
           if( kinv .ge. k0/4 ) then        
              id_exit = .true.
              if (scverbose) then
                call write_parallel('------- UW ShCu: Exit, kinv>k0/4')
              end if
              go to 333
           endif

!           kinv = min(kinv,7)
           ! From here, it must be '7 >= kinv >= 2'.        
 
       ! -------------------------------------------------------------------------- !
       ! Find PBL averaged tke ('tkeavg') and minimum 'thvl' ('thvlmin') in the PBL !
       ! In the current code, 'tkeavg' is obtained by averaging all interfacial TKE !
       ! within the PBL. However, in order to be conceptually consistent with   PBL !
       ! scheme, 'tkeavg' should be calculated by considering surface buoyancy flux.!
       ! If surface buoyancy flux is positive ( bflxs >0 ), surface interfacial TKE !
       ! should be included in calculating 'tkeavg', while if bflxs <= 0,   surface !
       ! interfacial TKE should not be included in calculating 'tkeavg'.   I should !
       ! modify the code when 'bflxs' is available as an input of cumulus scheme.   !
       ! 'thvlmin' is a minimum 'thvl' within PBL obtained by comparing top &  base !
       ! interface values of 'thvl' in each layers within the PBL.                  !
       ! -------------------------------------------------------------------------- !
       
         dpsum    = 0.
         thvlmin  = 1000.
         thvlavg  = 0.
         do k = 1,kinv ! max(kinv-1,1)    ! Here, 'k' is an interfacial layer index.  
            dpi = pifc0(k-1) - pifc0(k)
            dpsum  = dpsum  + dpi 
            tkeavg = tkeavg + dpi*tke(k)
            uavg   = uavg   + dpi*u0(k)
            vavg   = vavg   + dpi*v0(k)
            thvlavg = thvlavg + dpi*thvl0(k)
            if( k .ne. kinv ) thvlmin = min(thvlmin,min(thvl0bot(k),thvl0top(k)))
         end do
         tkeavg  = tkeavg/dpsum
         uavg    = uavg/dpsum
         vavg    = vavg/dpsum
         thvlavg = thvlavg/dpsum

        ! weighted average over lowest 20mb
!         dpsum = 0.
!         do k = 1,kinv
!             dpi = max(0.,(2e3+pmid0(k)-pifc0(0))/2e3)
!             qtavg  = qtavg  + dpi*qt0(k)
!             dpsum = dpsum + dpi
!         end do
!         qtavg   = qtavg/dpsum
 
       ! Interpolate qt to specified height
         k = 1
         do while (zmid0(k).lt.qtsrchgt)
           k = k+1
         end do
         if (k.gt.1) then
           qtavg = qt0(k-1)*(zmid0(k)-qtsrchgt) + qt0(k)*(qtsrchgt-zmid0(k-1))
           qtavg = qtavg / (zmid0(k)-zmid0(k-1))
         else
           qtavg = qt0(1)
         end if

       ! ------------------------------------------------------------------ !
       ! Find characteristics of cumulus source air: qtsrc,thlsrc,usrc,vsrc !
       ! Note that 'thlsrc' was concocted using 'thvlsrc' and 'qtsrc'.      !
       ! 'qtsrc' is defined as the lowest layer mid-point value;   'thlsrc' !
       ! is from 'qtsrc' and 'thvlmin=thvlsrc'; 'usrc' & 'vsrc' are defined !
       ! as the values just below the PBL top interface.                    !
       ! ------------------------------------------------------------------ !

         if (windsrcavg) then
            zrho = pifc0(0)/(287.04*(t0(1)*(1.+0.608*qv0(1))))
            buoyflx = (-shfx(i)/cp-0.608*t0(1)*evap(i))/zrho ! K m s-1
!            delzg = (zifc0(1)-zifc0(0))*g
            delzg = (50.0)*g   ! assume 50m surface scale
            wstar = max(0.,0.001-0.41*buoyflx*delzg/t0(1)) ! m3 s-3
            qpert_out(i) = 0.0
            tpert_out(i) = 0.0
            if (wstar > 0.001) then
              wstar = 1.0*wstar**.3333
              tpert_out(i) = thlsrc_fac*shfx(i)/(zrho*wstar*cp)  ! K
              qpert_out(i) = qtsrc_fac*evap(i)/(zrho*wstar)    ! kg kg-1
            end if
            qpert_out(i) = max(min(qpert_out(i),0.02*qt0(1)),0.)  ! limit to 1% of QT
            tpert_out(i) = 0.1+max(min(tpert_out(i),1.0),0.)          ! limit to 1K
            qtsrc   = qtavg + qpert_out(i)
!           qtsrc   = qt0(1) + qpert_out(i)
!           thvlsrc = thvlavg + tpert_out(i)*(1.0+zvir*qtsrc) !/exnmid0(1)
            thvlsrc = thvlmin + tpert_out(i)*(1.0+zvir*qtsrc) !/exnmid0(1)
            thlsrc  = thvlsrc / ( 1. + zvir * qtsrc )
            usrc  = uavg
            vsrc  = vavg
         else
            qtsrc   = qt0(1)
            thvlsrc = thvlmin
            thlsrc  = thvlsrc / ( 1. + zvir * qtsrc )
            usrc    = u0(kinv-1) + ssu0(kinv-1) * ( pifc0(kinv-1) - pmid0(kinv-1) )
            vsrc    = v0(kinv-1) + ssv0(kinv-1) * ( pifc0(kinv-1) - pmid0(kinv-1) )
         endif

         if (dotransport.eq.1) then
         do m = 1, ncnst
            trsrc(m) = tr0(1,m)
         enddo 
         endif

       ! ------------------------------------------------------------------ !
       ! Find LCL of the source air and a layer index containing LCL (klcl) !
       ! When the LCL is exactly at the interface, 'klcl' is a layer index  ! 
       ! having 'plcl' as the lower interface similar to the 'kinv' case.   !
       ! In the previous code, I assumed that if LCL is located within the  !
       ! lowest model layer ( 1 ) or the top model layer ( k0 ), then  no  !
       ! convective adjustment is performed and just exited.   However, in  !
       ! the revised code, I relaxed the first constraint and  even though  !
       ! LCL is at the lowest model layer, I allowed cumulus convection to  !
       ! be initiated. For this case, cumulus convection should be started  !
       ! from the PBL top height, as shown in the following code.           !
       ! When source air is already saturated even at the surface, klcl is  !
       ! set to 1.                                                          !
       ! ------------------------------------------------------------------ !

         if ( pifc0(0).lt.70000.or.pifc0(0).gt.115000. ) then
           if (scverbose) then
             print *,"pifc(0) outside valid range! pifc=",pifc0(0)
           end if
           id_exit = .true.
           go to 333
         end if
         if (qtsrc.gt.0.1.or.qtsrc.lt.1e-8 ) then
           print *,"qtsrc outside valid range! qtsrc=",qtsrc
           id_exit = .true.
           go to 333
         end if
         if (thlsrc.gt.400..or.thlsrc.lt.100. ) then
           print *,"thlsrc outside valid range! thlsrc=",thlsrc
           print *,"pifc0(0)=",pifc0(0)
           id_exit = .true.
           go to 333
        end if

         plcl = qsinvert(qtsrc,thlsrc,pifc0(0))
         do k = 0, k0
            if( pifc0(k) .lt. plcl ) then
                klcl = k
                go to 25
            endif
         end do
         klcl = k0
25       continue
         klcl = max(1,klcl)
     
         if( plcl .lt. 60000. ) then               
            exit_klclk0(i) = 1.
            id_exit = .true.
            if (scverbose) then
              call write_parallel('------- UW ShCu: exit, plcl<600mb')
            end if
            go to 333
         endif

        ! ------------------------------------------------------------- !
        ! Calculate environmental virtual potential temperature at LCL, !
        !'thv0lcl' which is solely used in the 'cin' calculation. Note  !
        ! that 'thv0lcl' is calculated first by calculating  'thl0lcl'  !
        ! and 'qt0lcl' at the LCL, and performing 'conden' afterward,   !
        ! in fully consistent with the other parts of the code.         !
        ! ------------------------------------------------------------- !

         thl0lcl = thl0(klcl) + ssthl0(klcl) * ( plcl - pmid0(klcl) )
         qt0lcl  = qt0(klcl)  + ssqt0(klcl)  * ( plcl - pmid0(klcl) )
         call conden(plcl,thl0lcl,qt0lcl,thj,qvj,qlj,qij,qse,id_check)
         if( id_check .eq. 1 ) then
            exit_conden(i) = 1.
            id_exit = .true.
            go to 333
         end if
         thv0lcl = thj * ( 1. + zvir * qvj - qlj - qij )


       ! ------------------------------------------------------------------------
       !
       ! Compute Convective Inhibition, 'cin' & 'cinlcl' [J/kg]=[m2/s2] TKE unit. !
       !
       !
       ! 'cin' (cinlcl) is computed from the PBL top interface to LFC (LCL) using ! 
       ! piecewisely reconstructed environmental profiles, assuming environmental !
       ! buoyancy profile within each layer ( or from LCL to upper interface in   !
       ! each layer ) is simply a linear profile. For the purpose of cin (cinlcl) !
       ! calculation, we simply assume that lateral entrainment does not occur in !
       ! updrafting cumulus plume, i.e., cumulus source air property is conserved.!
       ! Below explains some rules used in the calculations of cin (cinlcl).   In !
       ! general, both 'cin' and 'cinlcl' are calculated from a PBL top interface !
       ! to LCL and LFC, respectively :      !
       ! 1. If LCL is lower than the PBL height, cinlcl = 0 and cin is calculated !
       !    from PBL height to LFC.       !
       ! 2. If LCL is higher than PBL height,   'cinlcl' is calculated by summing !
       !    both positive and negative cloud buoyancy up to LCL using 'single_cin'!
       !    From the LCL to LFC, however, only negative cloud buoyancy is counted !
       !    to calculate final 'cin' upto LFC.       !
       ! 3. If either 'cin' or 'cinlcl' is negative, they are set to be zero.
       !
       ! In the below code, 'klfc' is the layer index containing 'LFC' similarto !
       ! 'kinv' and 'klcl'.
       !
       ! ------------------------------------------------------------------------!

        cin    = 0.
        cinlcl = 0.
        plfc   = 0.
        klfc   = k0



!        lts =  0.0
!        do k = 2,k0-1
!           if (pmid0(k).lt.70000.0) then
!             lts = t0(k-1)/exnmid0(k-1)  ! theta
!             exit
!           end if
!        end do
!        lts = lts - t0(1)/exnmid0(1)

        !------------------------------------------------------------------------!
        ! Case 1. LCL height is higher than PBL interface ( 'pLCL <=ps0(kinv-1)' ) !
        !------------------------------------------------------------------------!

        if( klcl .ge. kinv ) then

            do k = kinv, k0 - 1
               if( k .lt. klcl ) then
                   thvubot = thvlsrc
                   thvutop = thvlsrc  
                   cin     = cin + single_cin(pifc0(k-1),thv0bot(k),pifc0(k),thv0top(k),thvubot,thvutop)
               elseif( k .eq. klcl ) then
                   !----- Bottom to LCL
                   thvubot = thvlsrc
                   thvutop = thvlsrc
                   cin     = cin + single_cin(pifc0(k-1),thv0bot(k),plcl,thv0lcl,thvubot,thvutop)
                   if( cin .lt. 0. ) limit_cinlcl(i) = 1.
                   cinlcl  = max(cin,0.)
                   cin     = cinlcl
                   !----- LCL to Top
                   thvubot = thvlsrc
                   call conden(pifc0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check)
                   if( id_check .eq. 1 ) then
                       exit_conden(i) = 1.
                       id_exit = .true.
                       go to 333
                   end if
                   thvutop = thj * ( 1. + zvir*qvj - qlj - qij )
                   call getbuoy(plcl,thv0lcl,pifc0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                   if( plfc .gt. 0. ) then 
                       klfc = k 
                       go to 35
                   end if
               else
                   thvubot = thvutop
                   call conden(pifc0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check)
                   if( id_check .eq. 1 ) then
                       exit_conden(i) = 1.
                       id_exit = .true.
                       go to 333
                   end if
                   thvutop = thj * ( 1. + zvir*qvj - qlj - qij )
                   call getbuoy(pifc0(k-1),thv0bot(k),pifc0(k),thv0top(k),thvubot,thvutop,plfc,cin)
                   if( plfc .gt. 0. ) then 
                       klfc = k
                       go to 35
                   end if 
               endif
            end do        

       ! -----------------------------------------------------------------------!
       ! Case 2. LCL height is lower than PBL interface ( 'pLCL > ps0(kinv-1)') !
       ! -----------------------------------------------------------------------!

       else
          cinlcl = 0. 
          do k = kinv, k0 - 1
             call conden(pifc0(k-1),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check)
             if( id_check .eq. 1 ) then
                 exit_conden(i) = 1.
                 id_exit = .true.
                 go to 333
             end if
             thvubot = thj * ( 1. + zvir*qvj - qlj - qij )
             call conden(pifc0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check)
             if( id_check .eq. 1 ) then
                 exit_conden(i) = 1.
                 id_exit = .true.
                 go to 333
             end if
             thvutop = thj * ( 1. + zvir*qvj - qlj - qij )
             call getbuoy(pifc0(k-1),thv0bot(k),pifc0(k),thv0top(k),thvubot,thvutop,plfc,cin)
             if( plfc .gt. 0. ) then 
                 klfc = k
                 go to 35
             end if 
          end do
       endif  ! End of CIN case selection

 35    continue
       if( cin .lt. 0. ) limit_cin(i) = 1.
       cin = max(0.,cin)
!       cin = max(cin,0.04*(lts-18.))   ! kludge to reduce UW in StCu regions


       if( klfc .ge. k0 ) then
           klfc = k0
           if (scverbose) then
             call write_parallel('------ UWShCu: klfc >= k0')
           end if
           exit_klfck0(i) = 1.
           id_exit = .true.
           go to 333
       endif

       !----------------------------------------------------------------------!
       ! In order to calculate implicit 'cin' (or 'cinlcl'), save the initially !
       ! calculated 'cin' and 'cinlcl', and other related variables. These will !
       ! be restored after calculating implicit CIN.                            !
       !----------------------------------------------------------------------!

       if( iter .eq. 1 ) then 
           cin_i       = cin
           cinlcl_i    = cinlcl
           ke          = rbuoy / ( rkfre * tkeavg + epsvarw ) 
           kinv_o      = kinv     
           klcl_o      = klcl     
           klfc_o      = klfc    
           plcl_o      = plcl    
           plfc_o      = plfc     
           tkeavg_o    = tkeavg   
           thvlmin_o   = thvlmin
           qtsrc_o     = qtsrc    
           thvlsrc_o   = thvlsrc  
           thlsrc_o    = thlsrc
           usrc_o      = usrc     
           vsrc_o      = vsrc     
           thv0lcl_o   = thv0lcl  
           if (dotransport.eq.1) then
           do m = 1, ncnst
              trsrc_o(m) = trsrc(m)
           enddo
           end if
       endif


     ! Modification : If I impose w = max(0.1, w) up to the top interface of
     !                klfc, I should only use cinlfc.  That is, if I want to
     !                use cinlcl, I should not impose w = max(0.1, w).
     !                Using cinlcl is equivalent to treating only 'saturated'
     !                moist convection. Note that in this sense, I should keep
     !                the functionality of both cinlfc and cinlcl.
     !                However, the treatment of penetrative entrainment level becomes
     !                ambiguous if I choose 'cinlcl'. Thus, the best option is to use
     !                'cinlfc'.

       ! -------------------------------------------------------------------------- !
       ! Calculate implicit 'cin' by averaging initial and final cins.    Note that !
       ! implicit CIN is adopted only when cumulus convection stabilized the system,!
       ! i.e., only when 'del_CIN >0'. If 'del_CIN<=0', just use explicit CIN. Note !
       ! also that since 'cinlcl' is set to zero whenever LCL is below the PBL top, !
       ! (see above CIN calculation part), the use of 'implicit CIN=cinlcl'  is not !
       ! good. Thus, when using implicit CIN, always try to only use 'implicit CIN= !
       ! cin', not 'implicit CIN=cinlcl'. However, both 'CIN=cin' and 'CIN=cinlcl'  !
       ! are good when using explicit CIN.                                          !
       ! -------------------------------------------------------------------------- !

       if( iter .ne. 1 ) then

           cin_f = cin
           cinlcl_f = cinlcl
           if( use_CINcin ) then
               del_CIN = cin_f - cin_i
           else
               del_CIN = cinlcl_f - cinlcl_i
           endif

           if( del_CIN .gt. 0. ) then

               ! -------------------------------------------------------------- ! 
               ! Calculate implicit 'cin' and 'cinlcl'. Note that when we chose !
               ! to use 'implicit CIN = cin', choose 'cinlcl = cinlcl_i' below: !
               ! because iterative CIN only aims to obtain implicit CIN,  once  !
               ! we obtained 'implicit CIN=cin', it is good to use the original !
               ! profiles information for all the other variables after that.   !
               ! Note 'cinlcl' will be explicitly used in calculating  'wlcl' & !
               ! 'ufrclcl' after calculating 'winv' & 'ufrcinv'  at the PBL top !
               ! interface later, after calculating 'cbmf'.                     !
               ! -------------------------------------------------------------- !
         
               alpha = compute_alpha( del_CIN, ke ) 
               cin   = cin_i + alpha * del_CIN
               if( use_CINcin ) then
                   cinlcl = cinlcl_i                 
               else
                   cinlcl = cinlcl_i + alpha * del_cinlcl   
               endif

               ! ----------------------------------------------------------------- !
               ! Restore the original values from the previous 'iter_cin' step (1) !
               ! to compute correct tendencies for (n+1) time step by implicit CIN !
               ! ----------------------------------------------------------------- !

               kinv      = kinv_o     
               klcl      = klcl_o     
               klfc      = klfc_o    
               plcl      = plcl_o    
               plfc      = plfc_o     
               tkeavg    = tkeavg_o   
               thvlmin   = thvlmin_o
               qtsrc     = qtsrc_o    
               thvlsrc   = thvlsrc_o  
               thlsrc    = thlsrc_o
               usrc      = usrc_o     
               vsrc      = vsrc_o     
               thv0lcl   = thv0lcl_o  
               if (dotransport.eq.1) then
               do m = 1, ncnst
                  trsrc(m) = trsrc_o(m)
               enddo
               end if

               qv0(:k0)            = qv0_o(:k0)
               ql0(:k0)            = ql0_o(:k0)
               qi0(:k0)            = qi0_o(:k0)
               t0(:k0)             = t0_o(:k0)
               s0(:k0)             = s0_o(:k0)
               u0(:k0)             = u0_o(:k0)
               v0(:k0)             = v0_o(:k0)
               qt0(:k0)            = qt0_o(:k0)
               thl0(:k0)           = thl0_o(:k0)
               thvl0(:k0)          = thvl0_o(:k0)
               ssthl0(:k0)         = ssthl0_o(:k0)
               ssqt0(:k0)          = ssqt0_o(:k0)
               thv0bot(:k0)        = thv0bot_o(:k0)
               thv0top(:k0)        = thv0top_o(:k0)
               thvl0bot(:k0)       = thvl0bot_o(:k0)
               thvl0top(:k0)       = thvl0top_o(:k0)
               ssu0(:k0)           = ssu0_o(:k0) 
               ssv0(:k0)           = ssv0_o(:k0) 
               if (dotransport.eq.1) then
               do m = 1, ncnst
                  tr0(:k0,m)   = tr0_o(:k0,m)
                  sstr0(:k0,m) = sstr0_o(:k0,m)
               enddo
               end if

               ! ------------------------------------------------------ !
               ! Initialize all fluxes, tendencies, and other variables ! 
               ! in association with cumulus convection.                !
               ! ------------------------------------------------------ ! 

               umf(0:k0)          = 0.0
               dcm(:k0)           = 0.0
               emf(0:k0)          = 0.0
               slflx(0:k0)        = 0.0
               qtflx(0:k0)        = 0.0
               uflx(0:k0)         = 0.0
               vflx(0:k0)         = 0.0
               qvten(:k0)         = 0.0
               qlten(:k0)         = 0.0
               qiten(:k0)         = 0.0
               sten(:k0)          = 0.0
               uten(:k0)          = 0.0
               vten(:k0)          = 0.0
               qrten(:k0)         = 0.0
               qsten(:k0)         = 0.0
               dwten(:k0)         = 0.0
               diten(:k0)         = 0.0
               cufrc(:k0)         = 0.0
               qcu(:k0)           = 0.0
               qlu(:k0)           = 0.0
               qiu(:k0)           = 0.0
               fer(:k0)           = 0.0
               fdr(:k0)           = 0.0
               xco(:k0)           = 0.0
               qc(:k0)            = 0.0
!               qldet(:k0)         = 0.0
!               qidet(:k0)         = 0.0
               qlten_sub          = 0.0
               qiten_sub          = 0.0
               qc_l(:k0)          = 0.0
               qc_i(:k0)          = 0.0
               cbmf               = 0.0
               cnt                = k0
               cnb                = 0.0
               qtten(:k0)         = 0.0
               slten(:k0)         = 0.0
               ufrc(0:k0)         = 0.0

               thlu(0:k0)         = MAPL_UNDEF
               qtu(0:k0)          = MAPL_UNDEF
               uu(0:k0)           = MAPL_UNDEF
               vu(0:k0)           = MAPL_UNDEF
               wu(0:k0)           = MAPL_UNDEF
               thvu(0:k0)         = MAPL_UNDEF
               thlu_emf(0:k0)     = MAPL_UNDEF
               qtu_emf(0:k0)      = MAPL_UNDEF
               uu_emf(0:k0)       = MAPL_UNDEF
               vu_emf(0:k0)       = MAPL_UNDEF
             
               if (dotransport.eq.1) then
               do m = 1, ncnst
                  trflx(0:k0,m)   = 0.0
                  trten(:k0,m)    = 0.0
                  tru(0:k0,m)     = 0.0
                  tru_emf(0:k0,m) = 0.0
               enddo
               end if

               ! -------------------------------------------------- !
               ! Below are diagnostic output variables for detailed !
               ! analysis of cumulus scheme.                        !
               ! -------------------------------------------------- ! 

               ufrcinvbase        = 0.0
               ufrclcl            = 0.0
               winvbase           = 0.0
               wlcl               = 0.0
               emfkbup            = 0.0 
               cbmflimit          = 0.0
#ifdef UWDIAG
               excessu_arr(:k0)   = 0.0
               excess0_arr(:k0)   = 0.0
               xc_arr(:k0)        = 0.0
               aquad_arr(:k0)     = 0.0
               bquad_arr(:k0)     = 0.0
               cquad_arr(:k0)     = 0.0
               bogbot_arr(:k0)    = 0.0
               bogtop_arr(:k0)    = 0.0
#endif

          else ! When 'del_CIN < 0', use explicit CIN instead of implicit CIN.
           
               ! ----------------------------------------------------------- ! 
               ! Identifier showing whether explicit or implicit CIN is used !
               ! ----------------------------------------------------------- ! 

               ind_delcin(i) = 1.             
               if (scverbose) then
                 call write_parallel('------ UWShCu: del_CIN<0')
               end if

               ! --------------------------------------------------------- !
               ! Restore original output values of "iter_cin = 1" and exit !
               ! --------------------------------------------------------- !

               umf_out(i,0:k0)         = umf_s(0:k0)
               dcm_out(i,:k0)          = dcm_s(:k0)
               qvten_out(i,:k0)        = qvten_s(:k0)
               qlten_out(i,:k0)        = qlten_s(:k0)  
               qiten_out(i,:k0)        = qiten_s(:k0)
               sten_out(i,:k0)         = sten_s(:k0)
               uten_out(i,:k0)         = uten_s(:k0)  
               vten_out(i,:k0)         = vten_s(:k0)
               qrten_out(i,:k0)        = qrten_s(:k0)
               qsten_out(i,:k0)        = qsten_s(:k0)
               qldet_out(i,:k0)        = qldet_s(:k0)
               qidet_out(i,:k0)        = qidet_s(:k0)
               qlsub_out(i,:k0)        = qlsub_s(:k0)
               qisub_out(i,:k0)        = qisub_s(:k0)
               cush_inout(i)           = cush_s
               cufrc_out(i,:k0)        = cufrc_s(:k0)
               qtflx_out(i,0:k0)       = qtflx_s(0:k0)
               slflx_out(i,0:k0)       = slflx_s(0:k0)
               uflx_out(i,0:k0)        = uflx_s(0:k0)
               vflx_out(i,0:k0)        = vflx_s(0:k0)

#ifdef UWDIAG  
               qcu_out(i,:k0)          = qcu_s(:k0)    
               qlu_out(i,:k0)          = qlu_s(:k0)  
               qiu_out(i,:k0)          = qiu_s(:k0)  
               cbmf_out(i)             = cbmf_s
               qc_out(i,:k0)           = qc_s(:k0)  
               cnt_out(i)              = cnt_s
               cnb_out(i)              = cnb_s
               if (dotransport.eq.1) then
               do m = 1, ncnst
                  trten_out(i,:k0,m)   = trten_s(:k0,m)
               enddo  
               end if
#endif             

               ! ------------------------------------------------------------------------------ ! 
               ! Below are diagnostic output variables for detailed analysis of cumulus scheme. !
               ! The order of vertical index is reversed for this internal diagnostic output.   !
               ! ------------------------------------------------------------------------------ !   

               fer_out(i,1:k0)      = fer_s(:k0)  
               fdr_out(i,1:k0)      = fdr_s(:k0)  

#ifdef UWDIAG
               ufrcinvbase_out(i)       = ufrcinvbase_s
               ufrclcl_out(i)           = ufrclcl_s 
               winvbase_out(i)          = winvbase_s
               wlcl_out(i)              = wlcl_s
               plcl_out(i)              = plcl_s
               pinv_out(i)              = pinv_s    
               prel_out(i)              = prel_s    
               plfc_out(i)              = plfc_s    
               pbup_out(i)              = pbup_s
               ppen_out(i)              = ppen_s    
               qtsrc_out(i)             = qtsrc_s
               thlsrc_out(i)            = thlsrc_s
               thvlsrc_out(i)           = thvlsrc_s
               emfkbup_out(i)           = emfkbup_s
               cbmflimit_out(i)         = cbmflimit_s
               tkeavg_out(i)            = tkeavg_s
               zinv_out(i)              = zinv_s
               rcwp_out(i)              = rcwp_s
               rlwp_out(i)              = rlwp_s
               riwp_out(i)              = riwp_s

               xc_out(i,1:k0)           = xc_s(:k0)
               cinh_out(i)              = cin_s
               cinlclh_out(i)           = cinlcl_s

               wu_out(i,k0:0:-1)        = wu_s(0:k0)
               qtu_out(i,k0:0:-1)       = qtu_s(0:k0)
               thlu_out(i,k0:0:-1)      = thlu_s(0:k0)
               thvu_out(i,k0:0:-1)      = thvu_s(0:k0)
               uu_out(i,k0:0:-1)        = uu_s(0:k0)
               vu_out(i,k0:0:-1)        = vu_s(0:k0)
               qtu_emf_out(i,k0:0:-1)   = qtu_emf_s(0:k0)
               thlu_emf_out(i,k0:0:-1)  = thlu_emf_s(0:k0)
               uu_emf_out(i,k0:0:-1)    = uu_emf_s(0:k0)
               vu_emf_out(i,k0:0:-1)    = vu_emf_s(0:k0)
               uemf_out(i,k0:0:-1)      = uemf_s(0:k0)

               excessu_arr_out(i,k0:1:-1)  = excessu_arr_s(:k0)
               excess0_arr_out(i,k0:1:-1)  = excess0_arr_s(:k0)
               xc_arr_out(i,k0:1:-1)       = xc_arr_s(:k0)
               aquad_arr_out(i,k0:1:-1)    = aquad_arr_s(:k0)
               bquad_arr_out(i,k0:1:-1)    = bquad_arr_s(:k0)
               cquad_arr_out(i,k0:1:-1)    = cquad_arr_s(:k0)
               bogbot_arr_out(i,k0:1:-1)   = bogbot_arr_s(:k0)
               bogtop_arr_out(i,k0:1:-1)   = bogtop_arr_s(:k0)

               if (dotransport.eq.1) then
               do m = 1, ncnst
                  trflx_out(i,k0:0:-1,m)   = trflx_s(0:k0,m)  
                  tru_out(i,k0:0:-1,m)     = tru_s(0:k0,m)
                  tru_emf_out(i,k0:0:-1,m) = tru_emf_s(0:k0,m)
               enddo
               endif
#endif

              id_exit = .false.
              go to 333

           endif

         endif 

       ! ------------------------------------------------------------------ !
       ! Define a release level, 'prel' and release layer, 'krel'.          !
       ! 'prel' is the lowest level from which buoyancy sorting occurs, and !
       ! 'krel' is the layer index containing 'prel' in it, similar to  the !
       ! previous definitions of 'kinv', 'klcl', and 'klfc'.    In order to !
       ! ensure that only PBL scheme works within the PBL,  if LCL is below !
       ! PBL top height, then 'krel = kinv', while if LCL is above  PBL top !
       ! height, then 'krel = klcl'.   Note however that regardless of  the !
       ! definition of 'krel', cumulus convection induces fluxes within PBL !
       ! through 'fluxbelowinv'.  We can make cumulus convection start from !
       ! any level, even within the PBL by appropriately defining 'krel'  & !
       ! 'prel' here. Then it must be accompanied by appropriate definition !
       ! of source air properties, CIN, and re-setting of 'fluxbelowinv', & !
       ! many other stuffs.                                                 !
       ! Note that even when 'prel' is located above the PBL top height, we !
       ! still have cumulus convection between PBL top height and 'prel':   !
       ! we simply assume that no lateral mixing occurs in this range.      !
       ! ------------------------------------------------------------------ !

!         print *,'klcl=',klcl

         if( klcl .lt. kinv ) then
            krel    = kinv
            prel    = pifc0(krel-1)
            thv0rel = thv0bot(krel) 
         else
            krel    = klcl
            prel    = plcl 
            thv0rel = thv0lcl
         endif  

!         print *,'krel=',krel
       ! --------------------------------------------------------------------------- !
       ! Calculate cumulus base mass flux ('cbmf'), fractional area ('ufrcinv'), and !
       ! and mean vertical velocity (winv) of cumulus updraft at PBL top interface.  !
       ! Also, calculate updraft fractional area (ufrclcl) and vertical velocity  at !
       ! the LCL (wlcl). When LCL is below PBLH, cinlcl = 0 and 'ufrclcl = ufrcinv', !
       ! and 'wlcl = winv.                                                           !
       ! Only updrafts strong enough to overcome CIN can rise over PBL top interface.! 
       ! Thus,  in order to calculate cumulus mass flux at PBL top interface, 'cbmf',!
       ! we need to know 'CIN' ( the strength of potential energy barrier ) and      !
       ! 'sigmaw' ( a standard deviation of updraft vertical velocity at the PBL top !
       ! interface, a measure of turbulentce strength in the PBL ).   Naturally, the !
       ! ratio of these two variables, 'mu' - normalized CIN by TKE- is key variable !
       ! controlling 'cbmf'.  If 'mu' becomes large, only small fraction of updrafts !
       ! with very strong TKE can rise over the PBL - both 'cbmf' and 'ufrc' becomes !
       ! small, but 'winv' becomes large ( this can be easily understood by PDF of w !
       ! at PBL top ).  If 'mu' becomes small, lots of updraft can rise over the PBL !
       ! top - both 'cbmf' and 'ufrc' becomes large, but 'winv' becomes small. Thus, !
       ! all of the key variables associated with cumulus convection  at the PBL top !
       ! - 'cbmf', 'ufrc', 'winv' where 'cbmf = rho*ufrc*winv' - are a unique functi !
       ! ons of 'mu', normalized CIN. Although these are uniquely determined by 'mu',! 
       ! we usually impose two comstraints on 'cbmf' and 'ufrc': (1) because we will !
       ! simply assume that subsidence warming and drying of 'kinv-1' layer in assoc !
       ! iation with 'cbmf' at PBL top interface is confined only in 'kinv-1' layer, !
       ! cbmf must not be larger than the mass within the 'kinv-1' layer. Otherwise, !
       ! instability will occur due to the breaking of stability con. If we consider !
       ! semi-Lagrangian vertical advection scheme and explicitly consider the exten !
       ! t of vertical movement of each layer in association with cumulus mass flux, !
       ! we don't need to impose this constraint. However,  using a  semi-Lagrangian !
       ! scheme is a future research subject. Note that this constraint should be ap !
       ! plied for all interfaces above PBL top as well as PBL top interface.   As a !
       ! result, this 'cbmf' constraint impose a 'lower' limit on mu - 'mumin0'. (2) !
       ! in order for mass flux parameterization - rho*(w'a')= M*(a_c-a_e) - to   be !
       ! valid, cumulus updraft fractional area should be much smaller than 1.    In !
       ! current code, we impose 'rmaxfrac = 0.1 ~ 0.2'   through the whole vertical !
       ! layers where cumulus convection occurs. At the PBL top interface,  the same !
       ! constraint is made by imposing another lower 'lower' limit on mu, 'mumin1'. !
       ! After that, also limit 'ufrclcl' to be smaller than 'rmaxfrac' by 'mumin2'. !
       ! --------------------------------------------------------------------------- !
       
       ! --------------------------------------------------------------------------- !
       ! Calculate normalized CIN, 'mu' satisfying all the three constraints imposed !
       ! on 'cbmf'('mumin0'), 'ufrc' at the PBL top - 'ufrcinv' - ( by 'mumin1' from !
       ! a parameter sentence), and 'ufrc' at the LCL - 'ufrclcl' ( by 'mumin2').    !
       ! Note that 'cbmf' does not change between PBL top and LCL  because we assume !
       ! that buoyancy sorting does not occur when cumulus updraft is unsaturated.   !
       ! ---------------------------------------------------------------------------
 
         if( use_CINcin ) then       
            wcrit = sqrt( 2. * cin * rbuoy )      
         else
            wcrit = sqrt( 2. * cinlcl * rbuoy )   
         endif
         sigmaw = sqrt( rkfre * tkeavg + epsvarw )
         mu = wcrit/sigmaw/1.4142                  
         if( mu .ge. 3. ) then
            if (scverbose) then
              call write_parallel('mu >= 3')
            end if
            id_exit = .true.
            go to 333
         endif
         rho0inv = pifc0(kinv-1)/(r*thv0top(kinv-1)*exnifc0(kinv-1))
         cbmf = (rho0inv*sigmaw/2.5066)*exp(-mu**2)

         ! 1. 'cbmf' constraint
         cbmflimit = 0.9*dp0(kinv-1)/g/dt
         mumin0 = 0.
         if( cbmf .gt. cbmflimit ) mumin0 = sqrt(-log(2.5066*cbmflimit/rho0inv/sigmaw))
         ! 2. 'ufrcinv' constraint
         mu = max(max(mu,mumin0),mumin1)
         ! 3. 'ufrclcl' constraint      
         mulcl = sqrt(2.*cinlcl*rbuoy)/1.4142/sigmaw
         mulclstar = sqrt(max(0.,2.*(exp(-mu**2)/2.5066)**2*(1./erfc(mu)**2-0.25/rmaxfrac**2)))
         if( mulcl .gt. 1.e-8 .and. mulcl .gt. mulclstar ) then
            mumin2 = compute_mumin2(mulcl,rmaxfrac,mu)
            if( mu .gt. mumin2 ) then
                 call write_parallel('Critical error in mu calculation in UW_ShCu')
!                call endrun
            endif
            mu = max(mu,mumin2)
            if( mu .eq. mumin2 ) limit_ufrc(i) = 1.
         endif
         if( mu .eq. mumin0 ) limit_cbmf(i) = 1.
         if( mu .eq. mumin1 ) limit_ufrc(i) = 1.

       ! ------------------------------------------------------------------- !    
       ! Calculate final ['cbmf','ufrcinv','winv'] at the PBL top interface. !
       ! Note that final 'cbmf' here is obtained in such that 'ufrcinv' and  !
       ! 'ufrclcl' are smaller than ufrcmax with no instability.             !
       ! ------------------------------------------------------------------- !

         cbmf = 1.0*(rho0inv*sigmaw/2.5066)*exp(-mu**2)                       
         winv = sigmaw*(2./2.5066)*exp(-mu**2)/erfc(mu)
         ufrcinv = cbmf/winv/rho0inv


       ! ------------------------------------------------------------------- !
       ! Calculate ['ufrclcl','wlcl'] at the LCL. When LCL is below PBL top, !
       ! it automatically becomes 'ufrclcl = ufrcinv' & 'wlcl = winv', since !
       ! it was already set to 'cinlcl=0' if LCL is below PBL top interface. !
       ! Note 'cbmf' at the PBL top is the same as 'cbmf' at the LCL.  Note  !
       ! also that final 'cbmf' here is obtained in such that 'ufrcinv' and  !
       ! 'ufrclcl' are smaller than ufrcmax and there is no instability.     !
       ! By construction, it must be 'wlcl > 0' but for assurance, I checked !
       ! this again in the below block. If 'ufrclcl < 0.1%', just exit.      !
       ! ------------------------------------------------------------------- !

         wtw = winv * winv - 2. * cinlcl * rbuoy
         if( wtw .le. 0. ) then
            if (scverbose) then
              call write_parallel('wlcl < 0 at the LCL')
            end if
            exit_wtw(i) = 1.
            id_exit = .true.
            go to 333
         endif
         wlcl = sqrt(wtw)
         ufrclcl = cbmf/wlcl/rho0inv
         wrel = wlcl
         if( ufrclcl .le. 0.0001 ) then
            if (scverbose) then
              call write_parallel( 'ufrclcl <= 0.0001' ) 
            end if
            exit_ufrc(i) = 1.
            id_exit = .true.
            go to 333
         endif
         ufrc(krel-1) = ufrclcl

      ! ----------------------------------------------------------------------- !
       ! Below is just diagnostic output for detailed analysis of cumulus scheme !
       ! ----------------------------------------------------------------------- !

       ufrcinvbase        = ufrcinv
       winvbase           = winv
       umf(kinv-1:krel-1) = cbmf   
       wu(kinv-1:krel-1)  = winv

       ! -------------------------------------------------------------------------- ! 
       ! Define updraft properties at the level where buoyancy sorting starts to be !
       ! happening, i.e., by definition, at 'prel' level within the release layer.  !
       ! Because no lateral entrainment occurs upto 'prel', conservative scalars of ! 
       ! cumulus updraft at release level is same as those of source air.  However, ! 
       ! horizontal momentums of source air are modified by horizontal PGF forcings ! 
       ! from PBL top interface to 'prel'.  For this case, we should add additional !
       ! horizontal momentum from PBL top interface to 'prel' as will be done below !
       ! to 'usrc' and 'vsrc'. Note that below cumulus updraft properties - umf, wu,!
       ! thlu, qtu, thvu, uu, vu - are defined all interfaces not at the layer mid- !
       ! point. From the index notation of cumulus scheme, wu(k) is the cumulus up- !
       ! draft vertical velocity at the top interface of k layer.                   !
       ! Diabatic horizontal momentum forcing should be treated as a kind of 'body' !
       ! forcing without actual mass exchange between convective updraft and        !
       ! environment, but still taking horizontal momentum from the environment to  !
       ! the convective updrafts. Thus, diabatic convective momentum transport      !
       ! vertically redistributes environmental horizontal momentum.                !
       ! -------------------------------------------------------------------------- !

         emf(krel-1)  = 0.
         umf(krel-1)  = cbmf
         wu(krel-1)   = wrel
         thlu(krel-1) = thlsrc
         qtu(krel-1)  = qtsrc
         call conden(prel,thlsrc,qtsrc,thj,qvj,qlj,qij,qse,id_check)
         if( id_check .eq. 1 ) then
            exit_conden(i) = 1.
            id_exit = .true.
            if (scverbose) then
              call write_parallel('------- UW ShCu: exit, conden')
            end if
            go to 333
         endif
         thvu(krel-1) = thj * ( 1. + zvir*qvj - qlj - qij )       

         uplus = 0.
         vplus = 0.
         if( krel .eq. kinv ) then
            uplus = PGFc * ssu0(kinv) * ( prel - pifc0(kinv-1) )
            vplus = PGFc * ssv0(kinv) * ( prel - pifc0(kinv-1) )
         else
            do k = kinv, max(krel-1,kinv)
               uplus = uplus + PGFc * ssu0(k) * ( pifc0(k) - pifc0(k-1) )
               vplus = vplus + PGFc * ssv0(k) * ( pifc0(k) - pifc0(k-1) )
            end do
            uplus = uplus + PGFc * ssu0(krel) * ( prel - pifc0(krel-1) )
            vplus = vplus + PGFc * ssv0(krel) * ( prel - pifc0(krel-1) )
         end if
         uu(krel-1) = usrc + uplus
         vu(krel-1) = vsrc + vplus      

         if (dotransport.eq.1) then
         do m = 1, ncnst
            tru(krel-1,m)  = trsrc(m)
         enddo
         endif

       ! -------------------------------------------------------------------------- !
       ! Define environmental properties at the level where buoyancy sorting occurs !
       ! ('pe', normally, layer midpoint except in the 'krel' layer). In the 'krel' !
       ! layer where buoyancy sorting starts to occur, however, 'pe' is defined     !
       ! differently because LCL is regarded as lower interface for mixing purpose. !
       ! -------------------------------------------------------------------------- !

         pe      = 0.5 * ( prel + pifc0(krel) )
         qsat_pe = 0.5 * ( prel + pifc0(krel) )
         dpe     = prel - pifc0(krel)
         exne    = exnerfn(pe)
         thvebot = thv0rel
         thle    = thl0(krel) + ssthl0(krel) * ( pe - pmid0(krel) )
         qte     = qt0(krel)  + ssqt0(krel)  * ( pe - pmid0(krel) )
         ue      = u0(krel)   + ssu0(krel)   * ( pe - pmid0(krel) )
         ve      = v0(krel)   + ssv0(krel)   * ( pe - pmid0(krel) )
         if (dotransport.eq.1) then
         do m = 1, ncnst
            tre(m) = tr0(krel,m)  + sstr0(krel,m) * ( pe - pmid0(krel) )
         enddo
         end if

       !-------------------------! 
       ! Buoyancy-Sorting Mixing !
       !-------------------------!------------------------------------------------ !
       !                                                                           !
       !  In order to complete buoyancy-sorting mixing at layer mid-point, and so  ! 
       !  calculate 'updraft mass flux, updraft w velocity, conservative scalars'  !
       !  at the upper interface of each layer, we need following 3 information.   ! 
       !                                                                           !
       !  1. Pressure where mixing occurs ('pe'), and temperature at 'pe' which is !
       !     necessary to calculate various thermodynamic coefficients at pe. This !
       !     temperature is obtained by undiluted cumulus properties lifted to pe. ! 
       !  2. Undiluted updraft properties at pe - conservative scalar and vertical !
       !     velocity -which are assumed to be the same as the properties at lower !
       !     interface only for calculation of fractional lateral entrainment  and !
       !     detrainment rate ( fer(k) and fdr(k) [Pa-1] ), respectively.    Final !
       !     values of cumulus conservative scalars and w at the top interface are !
       !     calculated afterward after obtaining fer(k) & fdr(k).                 !
       !  3. Environmental properties at pe.                                       !
       ! ------------------------------------------------------------------------- !
       
       ! ------------------------------------------------------------------------ ! 
       ! Define cumulus scale height.                                             !
       ! Cumulus scale height is defined as the maximum height cumulus can reach. !
       ! In case of premitive code, cumulus scale height ('cush')  at the current !
       ! time step was assumed to be the same as 'cush' of previous time step.    !
       ! However, I directly calculated cush at each time step using an iterative !
       ! method. Note that within the cumulus scheme, 'cush' information is  used !
       ! only at two places during buoyancy-sorting process:                      !
       ! (1) Even negatively buoyancy mixtures with strong vertical velocity      !
       !     enough to rise up to 'rle*scaleh' (rle = 0.1) from pe are entrained  !
       !     into cumulus updraft,                                                !  
       ! (2) The amount of mass that is involved in buoyancy-sorting mixing       !
       !      process at pe is rei(k) = rkm/scaleh/rho*g [Pa-1]                   !
       ! In terms of (1), I think critical stopping distance might be replaced by !
       ! layer thickness. In future, we will use rei(k) = (0.5*rkm/z0(k)/rho/g).  !
       ! In the premitive code,  'scaleh' was largely responsible for the jumping !
       ! variation of precipitation amount.                                       !
       ! ------------------------------------------------------------------------ !

         scaleh = tscaleh
         if( tscaleh .le. 0.0 ) scaleh = 1000.


     ! Save time : Set iter_scaleh = 1. This will automatically use 'cush' from the previous time step
     !             at the first implicit iteration. At the second implicit iteration, it will use
     !             the updated 'cush' by the first implicit cin. So, this updating has an effect of
     !             doing one iteration for cush calculation, which is good. 
     !             So, only this setting of 'iter_scaleh = 1' is sufficient-enough to save computation time.
     ! OK

!         do iter_scaleh = 1, 3
         iter_scaleh = 1

       ! ---------------------------------------------------------------- !
       ! Initialization of 'kbup' and 'kpen'                              !
       ! ---------------------------------------------------------------- !
       ! 'kbup' is the top-most layer in which cloud buoyancy is positive !
       ! both at the top and bottom interface of the layer. 'kpen' is the !
       ! layer upto which cumulus panetrates ,i.e., cumulus w at the base !
       ! interface is positive, but becomes negative at the top interface.!
       ! Here, we initialize 'kbup' and 'kpen'. These initializations are !  
       ! not trivial but important, expecially   in calculating turbulent !
       ! fluxes without confliction among several physics as explained in !
       ! detail in the part of turbulent fluxes calculation later.   Note !
       ! that regardless of whether 'kbup' and 'kpen' are updated or  not !
       ! during updraft motion,  penetrative entrainments are dumped down !
       ! across the top interface of 'kbup' later.      More specifically,!
       ! penetrative entrainment heat and moisture fluxes are  calculated !
       ! from the top interface of 'kbup' layer  to the base interface of !
       ! 'kpen' layer. Because of this, initialization of 'kbup' & 'kpen' !
       ! influence the convection system when there are not updated.  The !  
       ! below initialization of 'kbup = krel' assures  that  penetrative !
       ! entrainment fluxes always occur at interfaces above the PBL  top !
       ! interfaces (i.e., only at interfaces k >=kinv ), which seems  to !
       ! be attractable considering that the most correct fluxes  at  the !
       ! PBL top interface can be ontained from the 'fluxbelowinv'  using !
       ! reconstructed PBL height.                                        ! 
       ! The 'kbup = krel'(after going through the whole buoyancy sorting !
       ! proces during updraft motion) implies that cumulus updraft  from !
       ! the PBL top interface can not reach to the LFC,so that 'kbup' is !
       ! not updated during upward. This means that cumulus updraft   did !
       ! not fully overcome the buoyancy barrier above just the PBL top.  !
       ! If 'kpen' is not updated either ( i.e., cumulus cannot rise over !
       ! the top interface of release layer),penetrative entrainment will !
       ! not happen at any interfaces.  If cumulus updraft can rise above !
       ! the release layer but cannot fully overcome the buoyancy barrier !
       ! just above PBL top interface, penetratve entrainment   occurs at !
       ! several above interfaces, including the top interface of release ! 
       ! layer. In the latter case, warming and drying tendencies will be !
       ! be initiated in 'krel' layer. Note current choice of 'kbup=krel' !
       ! is completely compatible with other flux physics without  double !
       ! or miss counting turbulent fluxes at any interface. However, the !
       ! alternative choice of 'kbup=krel-1' also has itw own advantage - !
       ! when cumulus updraft cannot overcome buoyancy barrier just above !
       ! PBL top, entrainment warming and drying are concentrated in  the !
       ! 'kinv-1' layer instead of 'kinv' layer for this case. This might !
       ! seems to be more dynamically reasonable, but I will choose the   !
       ! 'kbup = krel' choice since it is more compatible  with the other !
       ! parts of the code, expecially, when we chose ' use_emf=.false. ' !
       ! as explained in detail in turbulent flux calculation part.       !
       ! ---------------------------------------------------------------- ! 

         kbup    = krel
         kpen    = krel

       ! ------------------------------------------------------------ !
       ! Since 'wtw' is continuously updated during vertical motion,  !
       ! I need below initialization command within this 'iter_scaleh'!
       ! do loop. Similarily, I need initializations of environmental !
       ! properties at 'krel' layer as below.                         !
       ! ------------------------------------------------------------ !

         wtw     = wlcl * wlcl
         pe      = 0.5 * ( prel + pifc0(krel) )
         dpe     = prel - pifc0(krel)
         exne    = exnerfn(pe)
         thvebot = thv0rel
         thle    = thl0(krel) + ssthl0(krel) * ( pe - pmid0(krel) )
         qte     = qt0(krel)  + ssqt0(krel)  * ( pe - pmid0(krel) )
         ue      = u0(krel)   + ssu0(krel)   * ( pe - pmid0(krel) )
         ve      = v0(krel)   + ssv0(krel)   * ( pe - pmid0(krel) )
         if (dotransport.eq.1) then
         do m = 1, ncnst
            tre(m) = tr0(krel,m)  + sstr0(krel,m)  * ( pe - pmid0(krel) )
         enddo
         endif

       ! ----------------------------------------------------------------------- !
       ! Cumulus rises upward from 'prel' ( or base interface of  'krel' layer ) !
       ! until updraft vertical velocity becomes zero.                           !
       ! Buoyancy sorting is performed via two stages. (1) Using cumulus updraft !
       ! properties at the base interface of each layer,perform buoyancy sorting !
       ! at the layer mid-point, 'pe',  and update cumulus properties at the top !
       ! interface, and then  (2) by averaging updated cumulus properties at the !
       ! top interface and cumulus properties at the base interface,   calculate !
       ! cumulus updraft properties at pe that will be used  in buoyancy sorting !
       ! mixing - thlue, qtue and, wue.  Using this averaged properties, perform !
       ! buoyancy sorting again at pe, and re-calculate fer(k) and fdr(k). Using !
       ! this recalculated fer(k) and fdr(k),  finally calculate cumulus updraft !
       ! properties at the top interface - thlu, qtu, thvu, uu, vu. In the below,!
       ! 'iter_xc = 1' performs the first stage, while 'iter_xc= 2' performs the !
       ! second stage. We can increase the number of iterations, 'nter_xc'.as we !
       ! want, but a sample test indicated that about 3 - 5 iterations  produced !
       ! satisfactory converent solution. Finally, identify 'kbup' and 'kpen'.   !
       ! ----------------------------------------------------------------------- !

       do k = krel, k0 - 1 ! Here, 'k' is a layer index.

          km1 = k - 1

          thlue = thlu(km1)
          qtue  = qtu(km1)    
          wue   = wu(km1)
          wtwb  = wtw  

         do iter_xc = 1, niter_xc
          
            wtw = wu(km1) * wu(km1)

          ! ---------------------------------------------------------------- !
          ! Calculate environmental and cumulus saturation 'excess' at 'pe'. !
          ! Note that in order to calculate saturation excess, we should use ! 
          ! liquid water temperature instead of temperature  as the argument !
          ! of "qsat". But note normal argument of "qsat" is temperature.    !
          ! ---------------------------------------------------------------- !

            call conden(pe,thle,qte,thj,qvj,qlj,qij,qse,id_check)
            if( id_check .eq. 1 ) then
               exit_conden(i) = 1.
               id_exit = .true.
               if (scverbose) then
                 call write_parallel('------- UW ShCu: exit, conden')
               end if
               go to 333
            end if
            thv0j    = thj * ( 1. + zvir*qvj - qlj - qij )
            rhomid0j    = pe / ( r * thv0j * exne )
            qsat_arg = thle*exne
            qs =  GEOS_QSAT(qsat_arg,qsat_pe/100.)     
            excess0  = qte - qs

            call conden(pe,thlue,qtue,thj,qvj,qlj,qij,qse,id_check)
            if( id_check .eq. 1 ) then
               exit_conden(i) = 1.
               id_exit = .true.
               if (scverbose) then
                 call write_parallel('------- UW ShCu: exit, conden')
               end if
               go to 333
            end if
          ! ----------------------------------------------------------------- !
          ! Detrain excessive condensate larger than 'criqc' from the cumulus ! 
          ! updraft before performing buoyancy sorting. All I should to do is !
          ! to update 'thlue' &  'que' here. Below modification is completely !
          ! compatible with the other part of the code since 'thule' & 'qtue' !
          ! are used only for buoyancy sorting. I found that as long as I use !
          ! 'niter_xc >= 2',  detraining excessive condensate before buoyancy !
          ! sorting has negligible influence on the buoyancy sorting results. !   
          ! ----------------------------------------------------------------- !
            if( (qlj + qij) .gt. criqc ) then
               exql  = ( ( qlj + qij ) - criqc ) * qlj / ( qlj + qij )
               exqi  = ( ( qlj + qij ) - criqc ) * qij / ( qlj + qij )
               qtue  = qtue - exql - exqi
               thlue = thlue + (xlv/cp/exne)*exql + (xls/cp/exne)*exqi 
            endif
            call conden(pe,thlue,qtue,thj,qvj,qlj,qij,qse,id_check)
            if( id_check .eq. 1 ) then
               exit_conden(i) = 1.
               id_exit = .true.
               if (scverbose) then
                 call write_parallel('------- UW ShCu: exit, conden')
               end if
               go to 333
            end if
            thvj     = thj * ( 1. + zvir * qvj - qlj - qij )
            tj       = thj * exne ! This 'tj' is used for computing thermo. coeffs. below
            qsat_arg = thlue*exne
            qs =  GEOS_QSAT(qsat_arg,qsat_pe/100.)
            excessu  = qtue - qs

          ! ------------------------------------------------------------------- !
          ! Calculate critical mixing fraction, 'xc'. Mixture with mixing ratio !
          ! smaller than 'xc' will be entrained into cumulus updraft.  Both the !
          ! saturated updrafts with 'positive buoyancy' or 'negative buoyancy + ! 
          ! strong vertical velocity enough to rise certain threshold distance' !
          ! are kept into the updraft in the below program. If the core updraft !
          ! is unsaturated, we can set 'xc = 0' and let the cumulus  convection !
          ! still works or we may exit.                                         !
          ! Current below code does not entrain unsaturated mixture. However it !
          ! should be modified such that it also entrain unsaturated mixture.   !
          ! ------------------------------------------------------------------- !

          ! ----------------------------------------------------------------- !
          ! cridis : Critical stopping distance for buoyancy sorting purpose. !
          !          scaleh is only used here.                                !
          ! ----------------------------------------------------------------- !

           if (cridist_opt.eq.0) then
            cridis = rle*scaleh                 ! Original code
           else
            cridis = rle*(zifc0(k) - zifc0(k-1))  ! New code
           end if 

          ! ---------------- !
          ! Buoyancy Sorting !
          ! ---------------- !                   

          ! ----------------------------------------------------------------- !
          ! Case 1 : When both cumulus and env. are unsaturated or saturated. !
          ! ----------------------------------------------------------------- !
            xsat = 0.

            if( ( excessu .le. 0. .and. excess0 .le. 0. ) .or. ( excessu .ge. 0. .and. excess0 .ge. 0. ) ) then
                xc = min(1.,max(0.,1.-2.*rbuoy*g*cridis/wue**2.*(1.-thvj/thv0j)))
              ! Below 3 lines are diagnostic output not influencing
              ! numerical calculations.
                aquad = 0.
                bquad = 0.
                cquad = 0.
               if (excessu .gt. 0.) then
                  xsat = 1.
               else
                  xsat = 0.
               end if
            else
          ! -------------------------------------------------- !
          ! Case 2 : When either cumulus or env. is saturated. !
          ! -------------------------------------------------- !
              xsat    = excessu / ( excessu - excess0 );
              thlxsat = thlue + xsat * ( thle - thlue );
              qtxsat  = qtue  + xsat * ( qte - qtue );
              call conden(pe,thlxsat,qtxsat,thj,qvj,qlj,qij,qse,id_check)
              if( id_check .eq. 1 ) then
                  exit_conden(i) = 1.
                  id_exit = .true.
                  if (scverbose) then
                    call write_parallel('------- UW ShCu: exit, conden')
                  end if
                  go to 333
              end if
              thvxsat = thj * ( 1. + zvir * qvj - qlj - qij )               
              ! -------------------------------------------------- !
              ! kk=1 : Cumulus Segment, kk=2 : Environment Segment !
              ! -------------------------------------------------- ! 
              do kk = 1, 2 
                   if (xsat==1.) xsat=1.+1e-6
                   if( kk .eq. 1 ) then
                       thv_x0 = thvj;
                       thv_x1 = ( 1. - 1./xsat ) * thvj + ( 1./xsat ) * thvxsat;
                   else
                       thv_x1 = thv0j;
                       thv_x0 = ( xsat / ( xsat - 1. ) ) * thv0j + ( 1./( 1. - xsat ) ) * thvxsat;
                   endif
                   aquad =  wue**2;
                   bquad =  2.*rbuoy*g*cridis*(thv_x1 - thv_x0)/thv0j - 2.*wue**2;
                   cquad =  2.*rbuoy*g*cridis*(thv_x0 -  thv0j)/thv0j +    wue**2;
                   if( kk .eq. 1 ) then
                       if( ( bquad**2-4.*aquad*cquad ) .ge. 0. ) then
                             call roots(aquad,bquad,cquad,xs1,xs2,status)
                             x_cu = min(1.,max(0.,min(xsat,min(xs1,xs2))))
                       else
                             x_cu = xsat;
                       endif
                   else 
                       if( ( bquad**2-4.*aquad*cquad) .ge. 0. ) then
                             call roots(aquad,bquad,cquad,xs1,xs2,status)
                             x_en = min(1.,max(0.,max(xsat,min(xs1,xs2))))
                       else
                             x_en = 1.
                       endif
                   endif
              enddo
              if( x_cu .eq. xsat ) then
                  xc = max(x_cu, x_en)
              else
                  xc = x_cu
              endif
            endif

          ! ------------------------------------------------------------------------ !
          ! Compute fractional lateral entrainment & detrainment rate in each layers.!
          ! The unit of rei(k), fer(k), and fdr(k) is [Pa-1].  Alternative choice of !
          ! 'rei(k)' is also shown below, where coefficient 0.5 was from approximate !
          ! tuning against the BOMEX case.                                           !
          ! In order to prevent the onset of instability in association with cumulus !
          ! induced subsidence advection, cumulus mass flux at the top interface  in !
          ! any layer should be smaller than ( 90% of ) total mass within that layer.!
          ! I imposed limits on 'rei(k)' as below,  in such that stability condition ! 
          ! is always satisfied.                                                     !
          ! Below limiter of 'rei(k)' becomes negative for some cases, causing error.!
          ! So, for the time being, I came back to the original limiter.             !
          ! ------------------------------------------------------------------------ !
            ee2    = xc**2
            ud2    = 1. - 2.*xc + xc**2  ! (1-xc)**2
            if (min(scaleh,mixscale).ne.0.0) then
!              rei(k) = ( (rkm+max(0.,(zmid0(k)-detrhgt)/200.)) / min(scaleh,mixscale) / g / rhomid0j )   ! alternative
              rei(k) = ( (rkm+max(0.,(zmid0(k)-detrhgt)/200.)-max(0.,min(2.,(cnvtr(i))/2.5e-6))) / min(scaleh,mixscale) / g / rhomid0j )   ! alternative
            else
              rei(k) = ( 0.5 * rkm / min(1500.,zmid0(k)) / g /rhomid0j )       ! Jason-2_0 version
            end if

            if( xc .gt. 0.5 ) rei(k) = min(rei(k),0.9*log(dp0(k)/g/dt/umf(km1) + 1.)/dpe/(2.*xc-1.))
            fer(k) = rei(k) * ee2
            fdr(k) = rei(k) * ud2
            xco(k) = xc

          


          ! ------------------------------------------------------------------------------ !
          ! Iteration Start due to 'maxufrc' constraint [ ****************************** ] ! 
          ! ------------------------------------------------------------------------------ !

          ! -------------------------------------------------------------------------- !
          ! Calculate cumulus updraft mass flux and penetrative entrainment mass flux. !
          ! Note that  non-zero penetrative entrainment mass flux will be asigned only !
          ! to interfaces from the top interface of 'kbup' layer to the base interface !
          ! of 'kpen' layer as will be shown later.                                    !
          ! -------------------------------------------------------------------------- !

            umf(k) = umf(km1) * exp( dpe * ( fer(k) - fdr(k) ) )
            emf(k) = 0.
   
            dcm(k) = 0.5*(umf(k)+umf(km1))*rei(k)*dpe*min(1.,max(0.,xsat-xc))
!            dcm(k) = min(1.,max(0.,xsat-xc))

          ! --------------------------------------------------------- !
          ! Compute cumulus updraft properties at the top interface.  !
          ! Also use Tayler expansion in order to treat limiting case !
          ! --------------------------------------------------------- !

            if( fer(k)*dpe .lt. 1.e-4 ) then
              thlu(k) = thlu(km1) + ( thle + ssthl0(k) * dpe / 2. - thlu(km1) ) * fer(k) * dpe
              qtu(k)  =  qtu(km1) + ( qte  +  ssqt0(k) * dpe / 2. -  qtu(km1) ) * fer(k) * dpe
              uu(k)   =   uu(km1) + ( ue   +   ssu0(k) * dpe / 2. -   uu(km1) ) * fer(k) * dpe - PGFc * ssu0(k) * dpe
              vu(k)   =   vu(km1) + ( ve   +   ssv0(k) * dpe / 2. -   vu(km1) ) * fer(k) * dpe - PGFc * ssv0(k) * dpe
              if (dotransport.eq.1) then
              do m = 1, ncnst
                 tru(k,m)  =  tru(km1,m) + ( tre(m)  + sstr0(k,m) * dpe / 2.  -  tru(km1,m) ) * fer(k) * dpe
              enddo
              end if
            else
              thlu(k) = ( thle + ssthl0(k) / fer(k) - ssthl0(k) * dpe / 2. ) -          &
                        ( thle + ssthl0(k) * dpe / 2. - thlu(km1) + ssthl0(k) / fer(k) ) * exp(-fer(k) * dpe)
              qtu(k)  = ( qte  +  ssqt0(k) / fer(k) -  ssqt0(k) * dpe / 2. ) -          &  
                        ( qte  +  ssqt0(k) * dpe / 2. -  qtu(km1) +  ssqt0(k) / fer(k) ) * exp(-fer(k) * dpe)
              uu(k) =   ( ue + ( 1. - PGFc ) * ssu0(k) / fer(k) - ssu0(k) * dpe / 2. ) - &
                        ( ue +     ssu0(k) * dpe / 2. -   uu(km1) + ( 1. - PGFc ) * ssu0(k) / fer(k) ) * exp(-fer(k) * dpe)
              vu(k) =   ( ve + ( 1. - PGFc ) * ssv0(k) / fer(k) - ssv0(k) * dpe / 2. ) - &
                        ( ve +     ssv0(k) * dpe / 2. -   vu(km1) + ( 1. - PGFc ) * ssv0(k) / fer(k) ) * exp(-fer(k) * dpe)
              if (dotransport.eq.1) then
              do m = 1, ncnst
                 tru(k,m)  = ( tre(m)  + sstr0(k,m) / fer(k) - sstr0(k,m) * dpe / 2. ) - &  
                             ( tre(m)  + sstr0(k,m) * dpe / 2. - tru(km1,m) + sstr0(k,m) / fer(k) ) * exp(-fer(k) * dpe)
              enddo
              end if
            end if

          !------------------------------------------------------------------- !
          ! Expel some of cloud water and ice from cumulus  updraft at the top !
          ! interface.  Note that this is not 'detrainment' term  but a 'sink' !
          ! term of cumulus updraft qt ( or one part of 'source' term of  mean !
          ! environmental qt ). At this stage, as the most simplest choice, if !
          ! condensate amount within cumulus updraft is larger than a critical !
          ! value, 'criqc', expels the surplus condensate from cumulus updraft !
          ! to the environment. A certain fraction ( e.g., 'frc_sus' ) of this !
          ! expelled condesnate will be in a form that can be suspended in the !
          ! layer k where it was formed, while the other fraction, '1-frc_sus' ! 
          ! will be in a form of precipitatble (e.g.,can potentially fall down !
          ! across the base interface of layer k ). In turn we should describe !
          ! subsequent falling of precipitable condensate ('1-frc_sus') across !
          ! the base interface of the layer k, &  evaporation of precipitating !
          ! water in the below layer k-1 and associated evaporative cooling of !
          ! the later, k-1, and falling of 'non-evaporated precipitating water !
          ! ( which was initially formed in layer k ) and a newly-formed preci !
          ! pitable water in the layer, k-1', across the base interface of the !
          ! lower layer k-1.  Cloud microphysics should correctly describe all !
          ! of these process.  In a near future, I should significantly modify !
          ! this cloud microphysics, including precipitation-induced downdraft !
          ! also.                                                              !
          ! ------------------------------------------------------------------ !

            call conden(pifc0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check)
            if( id_check .eq. 1 ) then
              exit_conden(i) = 1.
              id_exit = .true.
              if (scverbose) then
                call write_parallel('------- UW ShCu: exit, conden')
              end if
              go to 333
            end if
            if( (qlj + qij) .gt. criqc ) then
               exql    = ( ( qlj + qij ) - criqc ) * qlj / ( qlj + qij )
               exqi    = ( ( qlj + qij ) - criqc ) * qij / ( qlj + qij )
               ! ---------------------------------------------------------------- !
               ! It is very important to re-update 'qtu' and 'thlu'  at the upper ! 
               ! interface after expelling condensate from cumulus updraft at the !
               ! top interface of the layer. As mentioned above, this is a 'sink' !
               ! of cumulus qt (or equivalently, a 'source' of environmentasl qt),!
               ! not a regular convective'detrainment'.                           !
               ! ---------------------------------------------------------------- !
               qtu(k)  = qtu(k) - exql - exqi
               thlu(k) = thlu(k) + (xlv/exnifc0(k)/cp)*exql + (xls/exnifc0(k)/cp)*exqi 
               ! ---------------------------------------------------------------- !
               ! Expelled cloud condensate into the environment from the updraft. ! 
               ! After all the calculation later, 'dwten' and 'diten' will have a !
               ! unit of [ kg/kg/s ], because it is a tendency of qt. Restoration !
               ! of 'dwten' and 'diten' to this correct unit through  multiplying !
               ! 'umf(k)*g/dp0(k)' will be performed later after finally updating !
               ! 'umf' using a 'rmaxfrac' constraint near the end of this updraft !
               ! buoyancy sorting loop.                                           !
               ! ---------------------------------------------------------------- !
               dwten(k) = exql   
               diten(k) = exqi
            else
               dwten(k) = 0.
               diten(k) = 0.
            endif

          ! ----------------------------------------------------------------- ! 
          ! Update 'thvu(k)' after detraining condensate from cumulus updraft.!
          ! ----------------------------------------------------------------- ! 
            call conden(pifc0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check)
            if( id_check .eq. 1 ) then
               exit_conden(i) = 1.
               id_exit = .true.
               if (scverbose) then
                 call write_parallel('------- UW ShCu: exit, conden')
               end if
               go to 333
            end if  
            thvu(k) = thj * ( 1. + zvir * qvj - qlj - qij )

          ! ----------------------------------------------------------- ! 
          ! Calculate updraft vertical velocity at the upper interface. !
          ! In order to calculate 'wtw' at the upper interface, we use  !
          ! 'wtw' at the lower interface. Note  'wtw'  is continuously  ! 
          ! updated as cumulus updraft rises.                           !
          ! ----------------------------------------------------------- !

            bogbot = rbuoy * ( thvu(km1) / thvebot  - 1. ) ! Cloud buoyancy at base interface
            bogtop = rbuoy * ( thvu(k) / thv0top(k) - 1. ) ! Cloud buoyancy at top  interface

            delbog = bogtop - bogbot
            drage  = fer(k) * ( 1. + rdrag )
            expfac = exp(-2.*drage*dpe)

            wtwb = wtw
            if( drage*dpe .gt. 1.e-3 ) then
              wtw = wtw*expfac + (delbog + (1.-expfac)*(bogbot + delbog/(-2.*drage*dpe)))/(rhomid0j*drage)
            else
              wtw = wtw + dpe * ( bogbot + bogtop ) / rhomid0j
            endif

        ! Force the plume rise at least to klfc of the undiluted plume.
        ! Because even the below is not complete, I decided not to include this.

        ! if( k .le. klfc ) then
        !     wtw = max( 1.e-2, wtw )
        ! endif 
         
          ! -------------------------------------------------------------- !
          ! Repeat 'iter_xc' iteration loop until 'iter_xc = niter_xc'.    !
          ! Also treat the case even when wtw < 0 at the 'kpen' interface. !
          ! -------------------------------------------------------------- !  
          
            if( wtw .gt. 0. ) then   
              thlue = 0.5 * ( thlu(km1) + thlu(k) )
              qtue  = 0.5 * ( qtu(km1)  +  qtu(k) )         
              wue   = 0.5 *   sqrt( max( wtwb + wtw, 0. ) )
            else
              go to 111
            endif 

         end do  ! iter_xc loop

     111 continue

          ! --------------------------------------------------------------------------- ! 
          ! Add the contribution of self-detrainment  to vertical variations of cumulus !
          ! updraft mass flux. The reason why we are trying to include self-detrainment !
          ! is as follows.  In current scheme,  vertical variation of updraft mass flux !
          ! is not fully consistent with the vertical variation of updraft vertical w.  !
          ! For example, within a given layer, let's assume that  cumulus w is positive !
          ! at the base interface, while negative at the top interface. This means that !
          ! cumulus updraft cannot reach to the top interface of the layer. However,    !
          ! cumulus updraft mass flux at the top interface is not zero according to the !
          ! vertical tendency equation of cumulus mass flux.   Ideally, cumulus updraft ! 
          ! mass flux at the top interface should be zero for this case. In order to    !
          ! assures that cumulus updraft mass flux goes to zero when cumulus updraft    ! 
          ! vertical velocity goes to zero, we are imposing self-detrainment term as    !
          ! below by considering layer-mean cloud buoyancy and cumulus updraft vertical !
          ! velocity square at the top interface. Use of auto-detrainment term will  be !
          ! determined by setting 'use_self_detrain=.true.' in the parameter sentence.  !
          ! --------------------------------------------------------------------------- !
     
            if( use_self_detrain ) then
              autodet = min( 0.5*g*(bogbot+bogtop)/(max(wtw,0.)+1.e-4), 0. ) 
              umf(k)  = umf(k) * exp( 0.637*(dpe/rhomid0j/g) * autodet )   
            end if      
            if( umf(k) .eq. 0. ) wtw = -1.

          ! -------------------------------------- !
          ! Below block is just a dignostic output !
          ! -------------------------------------- ! 

#ifdef UWDIAG
            excessu_arr(k) = excessu
            excess0_arr(k) = excess0
            xc_arr(k)      = xc
            aquad_arr(k)   = aquad
            bquad_arr(k)   = bquad
            cquad_arr(K)   = cquad
            bogbot_arr(k)  = bogbot
            bogtop_arr(k)  = bogtop
#endif

          ! ------------------------------------------------------------------- !
          ! 'kbup' is the upper most layer in which cloud buoyancy  is positive ! 
          ! both at the base and top interface.  'kpen' is the upper most layer !
          ! up to cumulus can reach. Usually, 'kpen' is located higher than the !
          ! 'kbup'. Note we initialized these by 'kbup = krel' & 'kpen = krel'. !
          ! As explained before, it is possible that only 'kpen' is updated,    !
          ! while 'kbup' keeps its initialization value. For this case, current !
          ! scheme will simply turns-off penetrative entrainment fluxes and use ! 
          ! normal buoyancy-sorting fluxes for 'kbup <= k <= kpen-1' interfaces,!
          ! in order to describe shallow continental cumulus convection.        !
          ! ------------------------------------------------------------------- !

          ! if( bogbot .gt. 0. .and. bogtop .gt. 0. ) then 
          ! if( bogtop .gt. 0. ) then          
            if( bogtop .gt. 0. .and. wtw .gt. 0. ) then 
              kbup = k
            end if

            if( wtw .le. 0. ) then
              kpen = k
              go to 45
            end if

            wu(k) = sqrt(wtw)
            if( wu(k) .gt. 100. ) then
              exit_wu(i) = 1.
              id_exit = .true.
              if (scverbose) then
                call write_parallel('------- UW ShCu: exited, wu>100')
              end if
              go to 333
            endif

          ! ---------------------------------------------------------------------------- !
          ! Iteration end due to 'rmaxfrac' constraint [ ***************************** ] ! 
          ! ---------------------------------------------------------------------------- !

          ! ---------------------------------------------------------------------- !
          ! Calculate updraft fractional area at the upper interface and set upper ! 
          ! limit to 'ufrc' by 'rmaxfrac'. In order to keep the consistency  among !
          ! ['ufrc','umf','wu (or wtw)'], if ufrc is limited by 'rmaxfrac', either !
          ! 'umf' or 'wu' should be changed. Although both 'umf' and 'wu (wtw)' at !
          ! the current upper interface are used for updating 'umf' & 'wu'  at the !
          ! next upper interface, 'umf' is a passive variable not influencing  the !
          ! buoyancy sorting process in contrast to 'wtw'. This is a reason why we !
          ! adjusted 'umf' instead of 'wtw'. In turn we updated 'fdr' here instead !
          ! of 'fer',  which guarantees  that all previously updated thermodynamic !
          ! variables at the upper interface before applying 'rmaxfrac' constraint !
          ! are already internally consistent,  even though 'ufrc'  is  limited by !
          ! 'rmaxfrac'. Thus, we don't need to go through interation loop again.If !
          ! If we update 'fer' however, we should go through above iteration loop. !
          ! ---------------------------------------------------------------------- !
            
          rhoifc0j  = pifc0(k) / ( r * 0.5 * ( thv0bot(k+1) + thv0top(k) )*exnifc0(k) )
          ufrc(k) = umf(k) / ( rhoifc0j * wu(k) )
          if( ufrc(k) .gt. rmaxfrac ) then
              limit_ufrc(i) = 1. 
              ufrc(k) = rmaxfrac
              umf(k)  = rmaxfrac * rhoifc0j * wu(k)
              fdr(k)  = fer(k) - log( umf(k) / umf(km1) ) / dpe
          endif

          ! ------------------------------------------------------------ !
          ! Update environmental properties for at the mid-point of next !
          ! upper layer for use in buoyancy sorting.                     !
          ! ------------------------------------------------------------ ! 

          pe      = pmid0(k+1)
          dpe     = dp0(k+1)
          exne    = exnmid0(k+1)
          thvebot = thv0bot(k+1)
          thle    = thl0(k+1)
          qte     = qt0(k+1)
          ue      = u0(k+1)
          ve      = v0(k+1) 
          if (dotransport.eq.1) then
          do m = 1, ncnst
             tre(m)  = tr0(k+1,m)
          enddo
          endif

         end do  ! k loop

       ! ------------------------------------------------------------------------------- !
       ! Up to this point, we finished all of buoyancy sorting processes from the 'krel' !
       ! layer to 'kpen' layer: at the top interface of individual layers, we calculated !
       ! updraft and penetrative mass fluxes [ umf(k) & emf(k) = 0 ], updraft fractional !
       ! area [ ufrc(k) ],  updraft vertical velocity [ wu(k) ],  updraft  thermodynamic !
       ! variables [thlu(k),qtu(k),uu(k),vu(k),thvu(k)]. In the layer,we also calculated !
       ! fractional entrainment-detrainment rate [ fer(k), fdr(k) ], and detrainment ten !
       ! dency of water and ice from cumulus updraft [ dwten(k), diten(k) ]. In addition,!
       ! we updated and identified 'krel' and 'kpen' layer index, if any.  In the 'kpen' !
       ! layer, we calculated everything mentioned above except the 'wu(k)' and 'ufrc(k)'!
       ! since a real value of updraft vertical velocity is not defined at the kpen  top !
       ! interface (note 'ufrc' at the top interface of layer is calculated from 'umf(k)'!
       ! and 'wu(k)'). As mentioned before, special treatment is required when 'kbup' is !
       ! not updated and so 'kbup = krel'.                                               !
       ! ------------------------------------------------------------------------------- !
       
       ! ------------------------------------------------------------------------------ !
       ! During the 'iter_scaleh' iteration loop, non-physical ( with non-zero values ) !
       ! values can remain in the variable arrays above (also 'including' in case of wu !
       ! and ufrc at the top interface) the 'kpen' layer. This can happen when the kpen !
       ! layer index identified from the 'iter_scaleh = 1' iteration loop is located at !
       ! above the kpen layer index identified from   'iter_scaleh = 3' iteration loop. !
       ! Thus, in the following calculations, we should only use the values in each     !
       ! variables only up to finally identified 'kpen' layer & 'kpen' interface except ! 
       ! 'wu' and 'ufrc' at the top interface of 'kpen' layer.    Note that in order to !
       ! prevent any problems due to these non-physical values, I re-initialized    the !
       ! values of [ umf(kpen:k0), emf(kpen:k0), dwten(kpen+1:k0), diten(kpen+1:k0),! 
       ! fer(kpen:k0), fdr(kpen+1:k0), ufrc(kpen:k0) ] to be zero after 'iter_scaleh'!
       ! do loop.                                                                       !
       ! ------------------------------------------------------------------------------ !
       
45      continue

       ! ------------------------------------------------------------------------------ !
       ! Calculate 'ppen( < 0 )', updarft penetrative distance from the lower interface !
       ! of 'kpen' layer. Note that bogbot & bogtop at the 'kpen' layer either when fer !
       ! is zero or non-zero was already calculated above.                              !
       ! It seems that below qudarature solving formula is valid only when bogbot < 0.  !
       ! Below solving equation is clearly wrong ! I should revise this !               !
       ! ------------------------------------------------------------------------------ ! 
            
         if( drage .eq. 0. ) then
           aquad =  ( bogtop - bogbot ) / ( pifc0(kpen) - pifc0(kpen-1) )
           bquad =  2. * bogbot
           cquad = -wu(kpen-1)**2 * rhomid0j
           call roots(aquad,bquad,cquad,xc1,xc2,status)
           if( status .eq. 0 ) then
               if( xc1 .le. 0. .and. xc2 .le. 0. ) then
                   ppen = max( xc1, xc2 )
                   ppen = min( 0.,max( -dp0(kpen), ppen ) )  
               elseif( xc1 .gt. 0. .and. xc2 .gt. 0. ) then
                   ppen = -dp0(kpen)
                   if (scverbose) then
                     call write_parallel('Warning : UW-Cumulus penetrates up to kpen interface')
                   end if
               else
                   ppen = min( xc1, xc2 )
                   ppen = min( 0.,max( -dp0(kpen), ppen ) )  
               endif
           else
               ppen = -dp0(kpen)
               if (scverbose) then
                 call write_parallel('Warning : UW-Cumulus penetrates up to kpen interface')
               end if
           endif       
         else 
           ppen = compute_ppen(wtwb,drage,bogbot,bogtop,rhomid0j,dp0(kpen))
         endif
         if( ppen .eq. -dp0(kpen) .or. ppen .eq. 0. ) limit_ppen(i) = 1.

       ! -------------------------------------------------------------------- !
       ! Re-calculate the amount of expelled condensate from cloud updraft    !
       ! at the cumulus top. This is necessary for refined calculations of    !
       ! bulk cloud microphysics at the cumulus top. Note that ppen < 0.   !
       ! In the below, I explicitly calculate 'thlu_top' & 'qtu_top' by       !
       ! using non-zero 'fer(kpen)'.                                          !    
       ! -------------------------------------------------------------------- !

         if( fer(kpen)*(-ppen) .lt. 1.e-4 ) then
           thlu_top = thlu(kpen-1) + ( thl0(kpen) + ssthl0(kpen) * (-ppen) / 2. - thlu(kpen-1) ) * fer(kpen) * (-ppen)
           qtu_top  =  qtu(kpen-1) + (  qt0(kpen) +  ssqt0(kpen) * (-ppen) / 2.  - qtu(kpen-1) ) * fer(kpen) * (-ppen)
         else
           thlu_top = ( thl0(kpen) + ssthl0(kpen) / fer(kpen) - ssthl0(kpen) * (-ppen) / 2. ) - &
                      ( thl0(kpen) + ssthl0(kpen) * (-ppen) / 2. - thlu(kpen-1) + ssthl0(kpen) / fer(kpen) ) &
                      * exp(-fer(kpen) * (-ppen))
           qtu_top  = ( qt0(kpen)  +  ssqt0(kpen) / fer(kpen) -  ssqt0(kpen) * (-ppen) / 2. ) - &  
                      ( qt0(kpen)  +  ssqt0(kpen) * (-ppen) / 2. -  qtu(kpen-1) +  ssqt0(kpen) / fer(kpen) ) &
                      * exp(-fer(kpen) * (-ppen))
         end if

         call conden(pifc0(kpen-1)+ppen,thlu_top,qtu_top,thj,qvj,qlj,qij,qse,id_check)
         if( id_check .eq. 1 ) then
           exit_conden(i) = 1.
           id_exit = .true.
           if (scverbose) then
             call write_parallel('------- UW ShCu: exit, conden')
           end if
           go to 333
         end if
         exntop = ((pifc0(kpen-1)+ppen)/p00)**rovcp
         if( (qlj + qij) .gt. criqc ) then
            dwten(kpen) = ( ( qlj + qij ) - criqc ) * qlj / ( qlj + qij )
            diten(kpen) = ( ( qlj + qij ) - criqc ) * qij / ( qlj + qij )
            qtu_top  = qtu_top - dwten(kpen) - diten(kpen)
            thlu_top = thlu_top + (xlv/cp/exntop)*dwten(kpen) + (xls/cp/exntop)*diten(kpen) 
         else
            dwten(kpen) = 0.
            diten(kpen) = 0.
         endif

       ! ----------------------------------------------------------------------- !
       ! Calculate cumulus scale height as the top height that cumulus can reach.!
       ! ----------------------------------------------------------------------- !
       
         rhoifc0j = pifc0(kpen-1)/(r*0.5*(thv0bot(kpen)+thv0top(kpen-1))*exnifc0(kpen-1))  
         cush   = zifc0(kpen-1) - ppen/rhoifc0j/g
         scaleh = cush 

!         end do  ! iter_scaleh loop

       ! -------------------------------------------------------------------- !   
       ! The 'forcedCu' is logical identifier saying whether cumulus updraft  !
       ! overcome the buoyancy barrier just above the PBL top. If it is true, !
       ! cumulus did not overcome the barrier -  this is a shallow convection !
       ! with negative cloud buoyancy, mimicking  shallow continental cumulus !
       ! convection. Depending on 'forcedCu' parameter, treatment of heat  &  !
       ! moisture fluxes at the entraining interfaces, 'kbup <= k < kpen - 1' !
       ! will be set up in a different ways, as will be shown later.          !
       ! -------------------------------------------------------------------- !
 
         if( kbup .eq. krel ) then 
           forcedCu = .true.
           limit_shcu(i) = 1.
         else
           forcedCu = .false.
           limit_shcu(i) = 0.
         endif  

       ! ------------------------------------------------------------------ !
       ! Filtering of unerasonable cumulus adjustment here.  This is a very !
       ! important process which should be done cautiously. Various ways of !
       ! filtering are possible depending on cases mainly using the indices !
       ! of key layers - 'klcl','kinv','krel','klfc','kbup','kpen'. At this !
       ! stage, the followings are all possible : 'kinv >= 2', 'klcl >= 1', !
       ! 'krel >= kinv', 'kbup >= krel', 'kpen >= krel'. I must design this !
       ! filtering very cautiously, in such that none of  realistic cumulus !
       ! convection is arbitrarily turned-off. Potentially, I might turn-off! 
       ! cumulus convection if layer-mean 'ql > 0' in the 'kinv-1' layer,in !
       ! order to suppress cumulus convection growing, based at the Sc top. ! 
       ! This is one of potential future modifications. Note that ppen < 0. !
       ! ------------------------------------------------------------------ !

         cldhgt = pifc0(kpen-1) + ppen
         if( forcedCu ) then
           if (scverbose) then
             call write_parallel( 'forcedCu - did not overcome initial buoyancy barrier')
           end if
           exit_cufilter(i) = 1.
           id_exit = .true.
           go to 333
         end if

       ! Limit 'additional shallow cumulus' for DYCOMS simulation.
       ! if( cldhgt.ge.88000. ) then
       !     id_exit = .true.
       !     go to 333
       ! end if
       
       ! ------------------------------------------------------------------------------ !
       ! Re-initializing some key variables above the 'kpen' layer in order to suppress !
       ! the influence of non-physical values above 'kpen', in association with the use !
       ! of 'iter_scaleh' loop. Note that umf, emf,  ufrc are defined at the interfaces !
       ! (0:k0), while 'dwten','diten', 'fer', 'fdr' are defined at layer mid-points.  !
       ! Initialization of 'fer' and 'fdr' is for correct writing purpose of diagnostic !
       ! output. Note that we set umf(kpen)=emf(kpen)=ufrc(kpen)=0, in consistent  with !
       ! wtw < 0  at the top interface of 'kpen' layer. However, we still have non-zero !
       ! expelled cloud condensate in the 'kpen' layer.                                 !
       ! ------------------------------------------------------------------------------ !

         umf(kpen:k0)     = 0.
         emf(kpen:k0)     = 0.
         ufrc(kpen:k0)    = 0.
         dwten(kpen+1:k0) = 0.
         diten(kpen+1:k0) = 0.
         fer(kpen+1:k0)   = 0.
         fdr(kpen+1:k0)   = 0.
         xco(kpen+1:k0)   = 0.

       ! ------------------------------------------------------------------------ !
       ! Calculate downward penetrative entrainment mass flux, 'emf(k) < 0',  and !
       ! thermodynamic properties of penetratively entrained airs at   entraining !
       ! interfaces. emf(k) is defined from the top interface of the  layer  kbup !
       ! to the bottom interface of the layer 'kpen'. Note even when  kbup = krel,!
       ! i.e.,even when 'kbup' was not updated in the above buoyancy  sorting  do !
       ! loop (i.e., 'kbup' remains as the initialization value),   below do loop !
       ! of penetrative entrainment flux can be performed without  any conceptual !
       ! or logical problems, because we have already computed all  the variables !
       ! necessary for performing below penetrative entrainment block.            !
       ! In the below 'do' loop, 'k' is an interface index at which non-zero 'emf'! 
       ! (penetrative entrainment mass flux) is calculated. Since cumulus updraft !
       ! is negatively buoyant in the layers between the top interface of 'kbup'  !
       ! layer (interface index, kbup) and the top interface of 'kpen' layer, the !
       ! fractional lateral entrainment, fer(k) within these layers will be close !
       ! to zero - so it is likely that only strong lateral detrainment occurs in !
       ! thses layers. Under this situation,we can easily calculate the amount of !
       ! detrainment cumulus air into these negatively buoyanct layers by  simply !
       ! comparing cumulus updraft mass fluxes between the base and top interface !
       ! of each layer: emf(k) = emf(k-1)*exp(-fdr(k)*dp0(k))                     !
       !                       ~ emf(k-1)*(1-rei(k)*dp0(k))                       !
       !                emf(k-1)-emf(k) ~ emf(k-1)*rei(k)*dp0(k)                  !
       ! Current code assumes that about 'rpen~10' times of these detrained  mass !
       ! are penetratively re-entrained down into the 'k-1' interface. And all of !
       ! these detrained masses are finally dumped down into the top interface of !
       ! 'kbup' layer. Thus, the amount of penetratively entrained air across the !
       ! top interface of 'kbup' layer with 'rpen~10' becomes too large.          !
       ! Note that this penetrative entrainment part can be completely turned-off !
       ! and we can simply use normal buoyancy-sorting involved turbulent  fluxes !
       ! by modifying 'penetrative entrainment fluxes' part below.                !
       ! ------------------------------------------------------------------------ !
       
       ! -----------------------------------------------------------------------!
       ! Calculate entrainment mass flux and conservative scalars of entraining !
       ! free air at interfaces of 'kbup <= k < kpen - 1'                       !
       ! ---------------------------------------------------------------------- !

         do k = 0, k0
           thlu_emf(k) = thlu(k)
           qtu_emf(k)  = qtu(k)
           uu_emf(k)   = uu(k)
           vu_emf(k)   = vu(k)
           if (dotransport.eq.1) then
           do m = 1, ncnst
             tru_emf(k,m)  = tru(k,m)
           enddo
           endif
         end do

         do k = kpen - 1, kbup, -1  ! Here, 'k' is an interface index at which
                                  ! penetrative entrainment fluxes are calculated. 
                                  
          rhoifc0j = pifc0(k) / ( r * 0.5 * ( thv0bot(k+1) + thv0top(k) ) * exnifc0(k) )

          if( k .eq. kpen - 1 ) then

             ! ------------------------------------------------------------------------ ! 
             ! Note that 'ppen' has already been calculated in the above 'iter_scaleh'  !
             ! loop assuming zero lateral entrainmentin the layer 'kpen'.               !
             ! ------------------------------------------------------------------------ !       
             
             ! -------------------------------------------------------------------- !
             ! Calculate returning mass flux, emf ( < 0 )                           !
             ! Current penetrative entrainment rate with 'rpen~10' is too large and !
             ! future refinement is necessary including the definition of 'thl','qt'! 
             ! of penetratively entrained air.  Penetratively entrained airs across !
             ! the 'kpen-1' interface is assumed to have the properties of the base !
             ! interface of 'kpen' layer. Note that 'emf ~ - umf/ufrc = - w * rho'. !
             ! Thus, below limit sets an upper limit of |emf| to be ~ 10cm/s, which !
             ! is very loose constraint. Here, I used more restricted constraint on !
             ! the limit of emf, assuming 'emf' cannot exceed a net mass within the !
             ! layer above the interface. Similar to the case of warming and drying !
             ! due to cumulus updraft induced compensating subsidence,  penetrative !
             ! entrainment induces compensating upwelling -     in order to prevent !  
             ! numerical instability in association with compensating upwelling, we !
             ! should similarily limit the amount of penetrative entrainment at the !
             ! interface by the amount of masses within the layer just above the    !
             ! penetratively entraining interface.                                  !
             ! -------------------------------------------------------------------- !
             
             if( ( umf(k)*ppen*rei(kpen)*rpen ) .lt. -0.1*rhoifc0j )         limit_emf(i) = 1.
             if( ( umf(k)*ppen*rei(kpen)*rpen ) .lt. -0.9*dp0(kpen)/g/dt ) limit_emf(i) = 1.             

             emf(k) = max( max( umf(k)*ppen*rei(kpen)*rpen, -0.1*rhoifc0j), -0.9*dp0(kpen)/g/dt)
             thlu_emf(k) = thl0(kpen) + ssthl0(kpen) * ( pifc0(k) - pmid0(kpen) )
             qtu_emf(k)  = qt0(kpen)  + ssqt0(kpen)  * ( pifc0(k) - pmid0(kpen) )
             uu_emf(k)   = u0(kpen)   + ssu0(kpen)   * ( pifc0(k) - pmid0(kpen) )     
             vu_emf(k)   = v0(kpen)   + ssv0(kpen)   * ( pifc0(k) - pmid0(kpen) )   
             if (dotransport.eq.1) then
             do m = 1, ncnst
                tru_emf(k,m)  = tr0(kpen,m)  + sstr0(kpen,m)  * ( pifc0(k) - pmid0(kpen) )
             enddo
             endif

          else ! if(k.lt.kpen-1). 
              
             ! --------------------------------------------------------------------------- !
             ! Note we are coming down from the higher interfaces to the lower interfaces. !
             ! Also note that 'emf < 0'. So, below operation is a summing not subtracting. !
             ! In order to ensure numerical stability, I imposed a modified correct limit  ! 
             ! of '-0.9*dp0(k+1)/g/dt' on emf(k).                                          !
             ! --------------------------------------------------------------------------- !

             if( use_cumpenent ) then  ! Original Cumulative Penetrative Entrainment

                 if( ( emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen ) .lt. -0.1*rhoifc0j )        limit_emf(i) = 1
                 if( ( emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen ) .lt. -0.9*dp0(k+1)/g/dt ) limit_emf(i) = 1         
                 emf(k) = max(max(emf(k+1)-umf(k)*dp0(k+1)*rei(k+1)*rpen, -0.1*rhoifc0j), -0.9*dp0(k+1)/g/dt )    
                 if( abs(emf(k)) .gt. abs(emf(k+1)) ) then
                     thlu_emf(k) = ( thlu_emf(k+1) * emf(k+1) + thl0(k+1) * ( emf(k) - emf(k+1) ) ) / emf(k)
                     qtu_emf(k)  = ( qtu_emf(k+1)  * emf(k+1) + qt0(k+1)  * ( emf(k) - emf(k+1) ) ) / emf(k)
                     uu_emf(k)   = ( uu_emf(k+1)   * emf(k+1) + u0(k+1)   * ( emf(k) - emf(k+1) ) ) / emf(k)
                     vu_emf(k)   = ( vu_emf(k+1)   * emf(k+1) + v0(k+1)   * ( emf(k) - emf(k+1) ) ) / emf(k)
                     if (dotransport.eq.1) then
                     do m = 1, ncnst
                        tru_emf(k,m)  = ( tru_emf(k+1,m)  * emf(k+1) + tr0(k+1,m)  * ( emf(k) - emf(k+1) ) ) / emf(k)
                     enddo
                     endif
                 else   
                     thlu_emf(k) = thl0(k+1)
                     qtu_emf(k)  =  qt0(k+1)
                     uu_emf(k)   =   u0(k+1)
                     vu_emf(k)   =   v0(k+1)
                     if (dotransport.eq.1) then
                     do m = 1, ncnst
                        tru_emf(k,m)  =  tr0(k+1,m)
                     enddo
                     endif
                 endif   
                     
             else ! Alternative Non-Cumulative Penetrative Entrainment

                 if( ( -umf(k)*dp0(k+1)*rei(k+1)*rpen ) .lt. -0.1*rhoifc0j )        limit_emf(i) = 1
                 if( ( -umf(k)*dp0(k+1)*rei(k+1)*rpen ) .lt. -0.9*dp0(k+1)/g/dt ) limit_emf(i) = 1         
                 emf(k) = max(max(-umf(k)*dp0(k+1)*rei(k+1)*rpen, -0.1*rhoifc0j), -0.9*dp0(k+1)/g/dt )    
                 thlu_emf(k) = thl0(k+1)
                 qtu_emf(k)  =  qt0(k+1)
                 uu_emf(k)   =   u0(k+1)
                 vu_emf(k)   =   v0(k+1)
                 if (dotransport.eq.1) then
                 do m = 1, ncnst
                    tru_emf(k,m)  =  tr0(k+1,m)
                 enddo
                 endif

             endif

          endif

          ! ---------------------------------------------------------------------------- !
          ! In this GCM modeling framework,  all what we should do is to calculate  heat !
          ! and moisture fluxes at the given geometrically-fixed height interfaces -  we !
          ! don't need to worry about movement of material height surface in association !
          ! with compensating subsidence or unwelling, in contrast to the bulk modeling. !
          ! In this geometrically fixed height coordinate system, heat and moisture flux !
          ! at the geometrically fixed height handle everything - a movement of material !
          ! surface is implicitly treated automatically. Note that in terms of turbulent !
          ! heat and moisture fluxes at model interfaces, both the cumulus updraft  mass !
          ! flux and penetratively entraining mass flux play the same role -both of them ! 
          ! warms and dries the 'kbup' layer, cools and moistens the 'kpen' layer,   and !
          ! cools and moistens any intervening layers between 'kbup' and 'kpen' layers.  !
          ! It is important to note these identical roles on turbulent heat and moisture !
          ! fluxes of 'umf' and 'emf'.                                                   !
          ! When 'kbup' is a stratocumulus-topped PBL top interface,  increase of 'rpen' !
          ! is likely to strongly diffuse stratocumulus top interface,  resulting in the !
          ! reduction of cloud fraction. In this sense, the 'kbup' interface has a  very !
          ! important meaning and role : across the 'kbup' interface, strong penetrative !
          ! entrainment occurs, thus any sharp gradient properties across that interface !
          ! are easily diffused through strong mass exchange. Thus, an initialization of ! 
          ! 'kbup' (and also 'kpen') should be done very cautiously as mentioned before. ! 
          ! In order to prevent this stron diffusion for the shallow cumulus convection  !
          ! based at the Sc top, it seems to be good to initialize 'kbup = krel', rather !
          ! that 'kbup = krel-1'.                                                        !
          ! ---------------------------------------------------------------------------- !
          
         end do  ! k loop

       !------------------------------------------------------------------ !
       !                                                                   ! 
       ! Compute turbulent heat, moisture, momentum flux at all interfaces !
       !                                                                   !
       !------------------------------------------------------------------ !
       ! It is very important to note that in calculating turbulent fluxes !
       ! below, we must not double count turbulent flux at any interefaces.!
       ! In the below, turbulent fluxes at the interfaces (interface index !
       ! k) are calculated by the following 4 blocks in consecutive order: !
       !                                                                   !
       ! (1) " 0 <= k <= kinv - 1 "  : PBL fluxes.                         !
       !     From 'fluxbelowinv' using reconstructed PBL height. Currently,!
       !     the reconstructed PBLs are independently calculated for  each !
       !     individual conservative scalar variables ( qt, thl, u, v ) in !
       !     each 'fluxbelowinv',  instead of being uniquely calculated by !
       !     using thvl. Turbulent flux at the surface is assumed to be 0. !
       ! (2) " kinv <= k <= krel - 1 " : Non-buoyancy sorting fluxes       !
       !     Assuming cumulus mass flux  and cumulus updraft thermodynamic !
       !     properties (except u, v which are modified by the PGFc during !
       !     upward motion) are conserved during a updraft motion from the !
       !     PBL top interface to the release level. If these layers don't !
       !     exist (e,g, when 'krel = kinv'), then  current routine do not !
       !     perform this routine automatically. So I don't need to modify !
       !     anything.                                                     ! 
       ! (3) " krel <= k <= kbup - 1 " : Buoyancy sorting fluxes           !
       !     From laterally entraining-detraining buoyancy sorting plumes. ! 
       ! (4) " kbup <= k < kpen-1 " : Penetrative entrainment fluxes       !
       !     From penetratively entraining plumes,                         !
       !                                                                   !
       ! In case of normal situation, turbulent interfaces  in each groups !
       ! are mutually independent of each other. Thus double flux counting !
       ! or ambiguous flux counting requiring the choice among the above 4 !
       ! groups do not occur normally. However, in case that cumulus plume !
       ! could not completely overcome the buoyancy barrier just above the !
       ! PBL top interface and so 'kbup = krel' (.forcedCu=.true.) ( here, !
       ! it can be either 'kpen = krel' as the initialization, or ' kpen > !
       ! krel' if cumulus updraft just penetrated over the top of  release !
       ! layer ). If this happens, we should be very careful in organizing !
       ! the sequence of the 4 calculation routines above -  note that the !
       ! routine located at the later has the higher priority.  Additional ! 
       ! feature I must consider is that when 'kbup = kinv - 1' (this is a !
       ! combined situation of 'kbup=krel-1' & 'krel = kinv' when I  chose !
       ! 'kbup=krel-1' instead of current choice of 'kbup=krel'), a strong !
       ! penetrative entrainment fluxes exists at the PBL top interface, & !
       ! all of these fluxes are concentrated (deposited) within the layer ! 
       ! just below PBL top interface (i.e., 'kinv-1' layer). On the other !
       ! hand, in case of 'fluxbelowinv', only the compensating subsidence !
       ! effect is concentrated in the 'kinv-1' layer and 'pure' turbulent !
       ! heat and moisture fluxes ( 'pure' means the fluxes not associated !
       ! with compensating subsidence) are linearly distributed throughout !
       ! the whole PBL. Thus different choice of the above flux groups can !
       ! produce very different results. Output variable should be written !
       ! consistently to the choice of computation sequences.              !
       ! When the case of 'kbup = krel(-1)' happens,another way to dealing !
       ! with this case is to simply ' exit ' the whole cumulus convection !
       ! calculation without performing any cumulus convection.     We can !
       ! choose this approach by specifying a condition in the  'Filtering !
       ! of unreasonable cumulus adjustment' just after 'iter_scaleh'. But !
       ! this seems not to be a good choice (although this choice was used !
       ! previous code ), since it might arbitrary damped-out  the shallow !
       ! cumulus convection over the continent land, where shallow cumulus ! 
       ! convection tends to be negatively buoyant.                        !
       ! ----------------------------------------------------------------- ! 

       ! --------------------------------------------------- !
       ! 1. PBL fluxes :  0 <= k <= kinv - 1                 !
       !    All the information necessary to reconstruct PBL ! 
       !    height are passed to 'fluxbelowinv'.             !
       ! --------------------------------------------------- !

         xsrc  = qtsrc
         xmean = qt0(kinv)
         xtop  = qt0(kinv+1) + ssqt0(kinv+1) * ( pifc0(kinv)   - pmid0(kinv+1) )
         xbot  = qt0(kinv-1) + ssqt0(kinv-1) * ( pifc0(kinv-1) - pmid0(kinv-1) )        
         call fluxbelowinv( cbmf, pifc0(0:k0), k0, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
         qtflx(0:kinv-1) = xflx(0:kinv-1)

         xsrc  = thlsrc
         xmean = thl0(kinv)
         xtop  = thl0(kinv+1) + ssthl0(kinv+1) * ( pifc0(kinv)   - pmid0(kinv+1) )
         xbot  = thl0(kinv-1) + ssthl0(kinv-1) * ( pifc0(kinv-1) - pmid0(kinv-1) )        
         call fluxbelowinv( cbmf, pifc0(0:k0), k0, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
         slflx(0:kinv-1) = cp * exnifc0(0:kinv-1) * xflx(0:kinv-1)

         xsrc  = usrc
         xmean = u0(kinv)
         xtop  = u0(kinv+1) + ssu0(kinv+1) * ( pifc0(kinv)   - pmid0(kinv+1) )
         xbot  = u0(kinv-1) + ssu0(kinv-1) * ( pifc0(kinv-1) - pmid0(kinv-1) )
         call fluxbelowinv( cbmf, pifc0(0:k0), k0, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
         uflx(0:kinv-1) = xflx(0:kinv-1)

         xsrc  = vsrc
         xmean = v0(kinv)
         xtop  = v0(kinv+1) + ssv0(kinv+1) * ( pifc0(kinv)   - pmid0(kinv+1) )
         xbot  = v0(kinv-1) + ssv0(kinv-1) * ( pifc0(kinv-1) - pmid0(kinv-1) )
         call fluxbelowinv( cbmf, pifc0(0:k0), k0, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
         vflx(0:kinv-1) = xflx(0:kinv-1)

         if (dotransport.eq.1) then
         do m = 1, ncnst
           xsrc  = trsrc(m)
           xmean = tr0(kinv,m)
           xtop  = tr0(kinv+1,m) + sstr0(kinv+1,m) * ( pifc0(kinv)   - pmid0(kinv+1) )
           xbot  = tr0(kinv-1,m) + sstr0(kinv-1,m) * ( pifc0(kinv-1) - pmid0(kinv-1) )        
           call fluxbelowinv( cbmf, pifc0(0:k0), k0, kinv, dt, xsrc, xmean, xtop, xbot, xflx )
           trflx(0:kinv-1,m) = xflx(0:kinv-1)
         enddo
         endif

       ! -------------------------------------------------------------- !
       ! 2. Non-buoyancy sorting fluxes : kinv <= k <= krel - 1         !
       !    Note that when 'krel = kinv', below block is never executed !
       !    as in a desirable, expected way ( but I must check  if this !
       !    is the case ). The non-buoyancy sorting fluxes are computed !
       !    only when 'krel > kinv'.                                    !
       ! -------------------------------------------------------------- !          

         uplus = 0.
         vplus = 0.
         do k = kinv, krel - 1
           kp1 = k + 1
           qtflx(k) = cbmf * ( qtsrc  - (  qt0(kp1) +  ssqt0(kp1) * ( pifc0(k) - pmid0(kp1) ) ) )          
           slflx(k) = cbmf * ( thlsrc - ( thl0(kp1) + ssthl0(kp1) * ( pifc0(k) - pmid0(kp1) ) ) ) * cp * exnifc0(k)
           uplus    = uplus + PGFc * ssu0(k) * ( pifc0(k) - pifc0(k-1) )
           vplus    = vplus + PGFc * ssv0(k) * ( pifc0(k) - pifc0(k-1) )
           uflx(k)  = cbmf * ( usrc + uplus -  (  u0(kp1)  +   ssu0(kp1) * ( pifc0(k) - pmid0(kp1) ) ) ) 
           vflx(k)  = cbmf * ( vsrc + vplus -  (  v0(kp1)  +   ssv0(kp1) * ( pifc0(k) - pmid0(kp1) ) ) )
           if (dotransport.eq.1) then
           do m = 1, ncnst
             trflx(k,m) = cbmf * ( trsrc(m)  - (  tr0(kp1,m) +  sstr0(kp1,m) * ( pifc0(k) - pmid0(kp1) ) ) )
           enddo          
           endif
         end do


       ! ------------------------------------------------------------------------ !
       ! 3. Buoyancy sorting fluxes : krel <= k <= kbup - 1                       !
       !    In case that 'kbup = krel - 1 ' ( or even in case 'kbup = krel' ),    ! 
       !    buoyancy sorting fluxes are not calculated, which is consistent,      !
       !    desirable feature.                                                    !  
       ! ------------------------------------------------------------------------ !

         do k = krel, kbup - 1      
           kp1 = k + 1
           slflx(k) = cp * exnifc0(k) * umf(k) * ( thlu(k) - ( thl0(kp1) + ssthl0(kp1) * ( pifc0(k) - pmid0(kp1) ) ) )
           qtflx(k) = umf(k) * ( qtu(k) - ( qt0(kp1) + ssqt0(kp1) * ( pifc0(k) - pmid0(kp1) ) ) )
           uflx(k)  = umf(k) * ( uu(k) - ( u0(kp1) + ssu0(kp1) * ( pifc0(k) - pmid0(kp1) ) ) )
           vflx(k)  = umf(k) * ( vu(k) - ( v0(kp1) + ssv0(kp1) * ( pifc0(k) - pmid0(kp1) ) ) )
           if (dotransport.eq.1) then
           do m = 1, ncnst
             trflx(k,m) = umf(k) * ( tru(k,m) - ( tr0(kp1,m) + sstr0(kp1,m) * ( pifc0(k) - pmid0(kp1) ) ) )
           enddo
           endif
         end do

       ! ------------------------------------------------------------------------- !
       ! 4. Penetrative entrainment fluxes : kbup <= k <= kpen - 1                 !
       !    The only confliction that can happen is when 'kbup = kinv-1'. For this !
       !    case, turbulent flux at kinv-1 is calculated  both from 'fluxbelowinv' !
       !    and here as penetrative entrainment fluxes.  Since penetrative flux is !
       !    calculated later, flux at 'kinv - 1 ' will be that of penetrative flux.!
       !    However, turbulent flux calculated at 'kinv - 1' from penetrative entr.!
       !    is less attractable,  since more reasonable turbulent flux at 'kinv-1' !
       !    should be obtained from 'fluxbelowinv', by considering  re-constructed ! 
       !    inversion base height. This conflicting problem can be solved if we can!
       !    initialize 'kbup = krel', instead of kbup = krel - 1. This choice seems!
       !    to be more reasonable since it is not conflicted with 'fluxbelowinv' in!
       !    calculating fluxes at 'kinv - 1' ( for this case, flux at 'kinv-1' is  !
       !    always from 'fluxbelowinv' ), and flux at 'krel-1' is calculated from  !
       !    the non-buoyancy sorting flux without being competed with penetrative  !
       !    entrainment fluxes. Even when we use normal cumulus flux instead of    !
       !    penetrative entrainment fluxes at 'kbup <= k <= kpen-1' interfaces,    !
       !    the initialization of kbup=krel perfectly works without any conceptual !
       !    confliction. Thus it seems to be much better to choose 'kbup = krel'   !
       !    initialization of 'kbup', which is current choice.                     !
       !    Note that below formula uses conventional updraft cumulus fluxes for   !
       !    shallow cumulus which did not overcome the first buoyancy barrier above!
       !    PBL top while uses penetrative entrainment fluxes for the other cases  !
       !    'kbup <= k <= kpen-1' interfaces. Depending on cases, however, I can   !
       !    selelct different choice.                                              !
       ! ------------------------------------------------------------------------------------------------------------------ !
       !   if( forcedCu ) then                                                                                              !
       !       slflx(k) = cp * exns0(k) * umf(k) * ( thlu(k) - ( thl0(kp1) + ssthl0(kp1) * ( ps0(k) - p0(kp1) ) ) )         !
       !       qtflx(k) =                 umf(k) * (  qtu(k) - (  qt0(kp1) +  ssqt0(kp1) * ( ps0(k) - p0(kp1) ) ) )         !
       !       uflx(k)  =                 umf(k) * (   uu(k) - (   u0(kp1) +   ssu0(kp1) * ( ps0(k) - p0(kp1) ) ) )         !
       !       vflx(k)  =                 umf(k) * (   vu(k) - (   v0(kp1) +   ssv0(kp1) * ( ps0(k) - p0(kp1) ) ) )         !
       !       do m = 1, ncnst                                                                                              !
       !          trflx(k,m) = umf(k) * ( tru(k,m) - ( tr0(kp1,m) + sstr0(kp1,m) * ( ps0(k) - p0(kp1) ) ) )                 !
       !       enddo                                                                                                        !
       !   else                                                                                                             !
       !       slflx(k) = cp * exns0(k) * emf(k) * ( thlu_emf(k) - ( thl0(k) + ssthl0(k) * ( ps0(k) - p0(k) ) ) )           !
       !       qtflx(k) =                 emf(k) * (  qtu_emf(k) - (  qt0(k) +  ssqt0(k) * ( ps0(k) - p0(k) ) ) )           !
       !       uflx(k)  =                 emf(k) * (   uu_emf(k) - (   u0(k) +   ssu0(k) * ( ps0(k) - p0(k) ) ) )           !
       !       vflx(k)  =                 emf(k) * (   vu_emf(k) - (   v0(k) +   ssv0(k) * ( ps0(k) - p0(k) ) ) )           !
       !       do m = 1, ncnst                                                                                              !
       !          trflx(k,m) = emf(k) * ( tru_emf(k,m) - ( tr0(k,m) + sstr0(k,m) * ( ps0(k) - p0(k) ) ) )                   !
       !       enddo                                                                                                        !
       !   endif                                                                                                            !
       !                                                                                                                    !
       !   if( use_uppenent ) then ! Combined Updraft + Penetrative Entrainment Flux                                        !
       !       slflx(k) = cp * exns0(k) * umf(k) * ( thlu(k)     - ( thl0(kp1) + ssthl0(kp1) * ( ps0(k) - p0(kp1) ) ) ) + & !
       !                  cp * exns0(k) * emf(k) * ( thlu_emf(k) - (   thl0(k) +   ssthl0(k) * ( ps0(k) - p0(k) ) ) )       !
       !       qtflx(k) =                 umf(k) * (  qtu(k)     - (  qt0(kp1) +  ssqt0(kp1) * ( ps0(k) - p0(kp1) ) ) ) + & !
       !                                  emf(k) * (  qtu_emf(k) - (    qt0(k) +    ssqt0(k) * ( ps0(k) - p0(k) ) ) )       !
       !       uflx(k)  =                 umf(k) * (   uu(k)     - (   u0(kp1) +   ssu0(kp1) * ( ps0(k) - p0(kp1) ) ) ) + & !
       !                                  emf(k) * (   uu_emf(k) - (     u0(k) +     ssu0(k) * ( ps0(k) - p0(k) ) ) )       !
       !       vflx(k)  =                 umf(k) * (   vu(k)     - (   v0(kp1) +   ssv0(kp1) * ( ps0(k) - p0(kp1) ) ) ) + & !
       !                                  emf(k) * (   vu_emf(k) - (     v0(k) +     ssv0(k) * ( ps0(k) - p0(k) ) ) )       !
       !       do m = 1, ncnst                                                                                              !
       !          trflx(k,m) = umf(k) * ( tru(k,m) - ( tr0(kp1,m) + sstr0(kp1,m) * ( ps0(k) - p0(kp1) ) ) ) + &             ! 
       !                       emf(k) * ( tru_emf(k,m) - ( tr0(k,m) + sstr0(k,m) * ( ps0(k) - p0(k) ) ) )                   ! 
       !       enddo                                                                                                        !
       ! ------------------------------------------------------------------------------------------------------------------ !

         do k = kbup, kpen - 1      
           kp1 = k + 1
           slflx(k) = cp * exnifc0(k) * emf(k) * ( thlu_emf(k) - ( thl0(k) + ssthl0(k) * ( pifc0(k) - pmid0(k) ) ) )
           qtflx(k) =                 emf(k) * (  qtu_emf(k) - (  qt0(k) +  ssqt0(k) * ( pifc0(k) - pmid0(k) ) ) ) 
           uflx(k)  =                 emf(k) * (   uu_emf(k) - (   u0(k) +   ssu0(k) * ( pifc0(k) - pmid0(k) ) ) ) 
           vflx(k)  =                 emf(k) * (   vu_emf(k) - (   v0(k) +   ssv0(k) * ( pifc0(k) - pmid0(k) ) ) )
           if (dotransport.eq.1) then
           do m = 1, ncnst
             trflx(k,m) = emf(k) * ( tru_emf(k,m) - ( tr0(k,m) + sstr0(k,m) * ( pifc0(k) - pmid0(k) ) ) ) 
           enddo
           endif
         end do

       ! ------------------------------------------- !
       ! Turn-off cumulus momentum flux as an option !
       ! ------------------------------------------- !

         if( .not. use_momenflx ) then
           uflx(0:k0) = 0.
           vflx(0:k0) = 0.
         endif  

       ! -------------------------------------------------------- !
       ! Condensate tendency by compensating subsidence/upwelling !
       ! -------------------------------------------------------- !
       
         uemf(0:k0)         = 0.
         do k = 0, kinv - 2  ! Assume linear updraft mass flux within the PBL.
           uemf(k) = cbmf * ( pifc0(0) - pifc0(k) ) / ( pifc0(0) - pifc0(kinv-1) ) 
         end do
         uemf(kinv-1:krel-1) = cbmf
         uemf(krel:kbup-1)   = umf(krel:kbup-1)
         uemf(kbup:kpen-1)   = emf(kbup:kpen-1) ! Only use penetrative entrainment flux consistently.

         comsub(1:k0) = 0.
         do k = 1, kpen
           comsub(k)  = 0.5 * ( uemf(k) + uemf(k-1) )  ! comsub defined on interfaces 
         end do     

         do k = 1, kpen
           if( comsub(k) .ge. 0. ) then
              if( k .eq. k0 ) then
                  thlten_sub = 0.
                  qtten_sub  = 0.
                  qlten_sub  = 0.
                  qiten_sub  = 0.
                  nlten_sub  = 0.
                  niten_sub  = 0.
              else
                  thlten_sub = g * comsub(k) * ( thl0(k+1) - thl0(k) ) / ( pmid0(k) - pmid0(k+1) )
                  qtten_sub  = g * comsub(k) * (  qt0(k+1) -  qt0(k) ) / ( pmid0(k) - pmid0(k+1) )
                  qlten_sub  = g * comsub(k) * (  ql0(k+1) -  ql0(k) ) / ( pmid0(k) - pmid0(k+1) )
                  qiten_sub  = g * comsub(k) * (  qi0(k+1) -  qi0(k) ) / ( pmid0(k) - pmid0(k+1) )
!                  nlten_sub  = g * comsub(k) * (  tr0(k+1,ixnumliq) -  tr0(k,ixnumliq) ) / ( pmid0(k) - pmid0(k+1) )
!                  niten_sub  = g * comsub(k) * (  tr0(k+1,ixnumice) -  tr0(k,ixnumice) ) / ( pmid0(k) - pmid0(k+1) )
              endif
           else
              if( k .eq. 1 ) then
                  thlten_sub = 0.
                  qtten_sub  = 0.
                  qlten_sub  = 0.
                  qiten_sub  = 0.
                  nlten_sub  = 0.
                  niten_sub  = 0.
              else
                  thlten_sub = g * comsub(k) * ( thl0(k) - thl0(k-1) ) / ( pmid0(k-1) - pmid0(k) )
                  qtten_sub  = g * comsub(k) * (  qt0(k) -  qt0(k-1) ) / ( pmid0(k-1) - pmid0(k) )
                  qlten_sub  = g * comsub(k) * (  ql0(k) -  ql0(k-1) ) / ( pmid0(k-1) - pmid0(k) )
                  qiten_sub  = g * comsub(k) * (  qi0(k) -  qi0(k-1) ) / ( pmid0(k-1) - pmid0(k) )
!                  nlten_sub  = g * comsub(k) * (  tr0(k,ixnumliq) -  tr0(k-1,ixnumliq) ) / ( pmid0(k-1) - pmid0(k) )
!                  niten_sub  = g * comsub(k) * (  tr0(k,ixnumice) -  tr0(k-1,ixnumice) ) / ( pmid0(k-1) - pmid0(k) )
              endif
           endif
           thl_prog = thl0(k) + thlten_sub * dt
           qt_prog  = max( qt0(k) + qtten_sub * dt, 1.e-12 )
           call conden(pmid0(k),thl_prog,qt_prog,thj,qvj,qlj,qij,qse,id_check)
           if( id_check .eq. 1 ) then
              id_exit = .true.
              if (scverbose) then
                call write_parallel('------- UW ShCu: exit, conden L3302')
              end if
              go to 333
           endif
         ! qlten_sink(k) = ( qlj - ql0(k) ) / dt
         ! qiten_sink(k) = ( qij - qi0(k) ) / dt
           qlten_sink(k) = max( qlten_sub, - ql0(k) / dt ) ! For consistency with prognostic macrophysics scheme
           qiten_sink(k) = max( qiten_sub, - qi0(k) / dt ) ! For consistency with prognostic macrophysics scheme
!           nlten_sink(k) = max( nlten_sub, - tr0(k,ixnumliq) / dt ) 
!           niten_sink(k) = max( niten_sub, - tr0(k,ixnumice) / dt )
         end do

       ! --------------------------------------------- !
       !                                               !
       ! Calculate convective tendencies at each layer ! 
       !                                               !
       ! --------------------------------------------- !
       
       ! ----------------- !
       ! Momentum tendency !
       ! ----------------- !
       
         do k = 1, kpen
           km1 = k - 1 
           uten(k) = ( uflx(km1) - uflx(k) ) * g / dp0(k)
           vten(k) = ( vflx(km1) - vflx(k) ) * g / dp0(k) 
           uf(k)   = u0(k) + uten(k) * dt
           vf(k)   = v0(k) + vten(k) * dt
         ! do m = 1, ncnst
         !    trten(k,m) = ( trflx(km1,m) - trflx(k,m) ) * g / dp0(k)
         !  ! Limit trten(k,m) such that negative value is not developed.
         !  ! This limitation does not conserve grid-mean tracers and future
         !  ! refinement is required for tracer-conserving treatment.
         !    trten(k,m) = max(trten(k,m),-tr0(k,m)/dt)              
         ! enddo
         end do  

       ! ----------------------------------------------------------------- !
       ! Tendencies of thermodynamic variables.                            ! 
       ! This part requires a careful treatment of bulk cloud microphysics.!
       ! Relocations of 'precipitable condensates' either into the surface ! 
       ! or into the tendency of 'krel' layer will be performed just after !
       ! finishing the below 'do-loop'.                                    !        
       ! ----------------------------------------------------------------- !
       

       do k = 1, kpen

          km1 = k - 1

          ! ------------------------------------------------------------------------------ !
          ! Compute 'slten', 'qtten', 'qvten', 'qlten', 'qiten', and 'sten'                !
          !                                                                                !
          ! Key assumptions made in this 'cumulus scheme' are :                            !
          ! 1. Cumulus updraft expels condensate into the environment at the top interface !
          !    of each layer. Note that in addition to this expel process ('source' term), !
          !    cumulus updraft can modify layer mean condensate through normal detrainment !
          !    forcing or compensating subsidence.                                         !
          ! 2. Expelled water can be either 'sustaining' or 'precipitating' condensate. By !
          !    definition, 'suataining condensate' will remain in the layer where it was   !
          !    formed, while 'precipitating condensate' will fall across the base of the   !
          !    layer where it was formed.                                                  !
          ! 3. All precipitating condensates are assumed to fall into the release layer or !
          !    ground as soon as it was formed without being evaporated during the falling !
          !    process down to the desinated layer ( either release layer of surface ).    !
          ! ------------------------------------------------------------------------------ !

          ! ------------------------------------------------------------------------- !     
          ! 'dwten(k)','diten(k)' : Production rate of condensate  within the layer k !
          !      [ kg/kg/s ]        by the expels of condensate from cumulus updraft. !
          ! It is important to note that in terms of moisture tendency equation, this !
          ! is a 'source' term of enviromental 'qt'.  More importantly,  these source !
          ! are already counted in the turbulent heat and moisture fluxes we computed !
          ! until now, assuming all the expelled condensate remain in the layer where ! 
          ! it was formed. Thus, in calculation of 'qtten' and 'slten' below, we MUST !
          ! NOT add or subtract these terms explicitly in order not to double or miss !
          ! count, unless some expelled condensates fall down out of the layer.  Note !
          ! this falling-down process ( i.e., precipitation process ) and  associated !
          ! 'qtten' and 'slten' and production of surface precipitation flux  will be !
          ! treated later in 'zm_conv_evap' in 'convect_shallow_tend' subroutine.     ! 
          ! In below, we are converting expelled cloud condensate into correct unit.  !
          ! I found that below use of '0.5 * (umf(k-1) + umf(k))' causes conservation !
          ! errors at some columns in global simulation. So, I returned to originals. !
          ! This will cause no precipitation flux at 'kpen' layer since umf(kpen)=0.  !
          ! ------------------------------------------------------------------------- !

            dwten(k) = dwten(k) * 0.5 * ( umf(k-1) + umf(k) ) * g / dp0(k) ! [ kg/kg/s ]
            diten(k) = diten(k) * 0.5 * ( umf(k-1) + umf(k) ) * g / dp0(k) ! [ kg/kg/s ] 

          ! dwten(k) = dwten(k) * umf(k) * g / dp0(k) ! [ kg/kg/s ]
          ! diten(k) = diten(k) * umf(k) * g / dp0(k) ! [ kg/kg/s ]

          ! --------------------------------------------------------------------------- !
          ! 'qrten(k)','qsten(k)' : Production rate of rain and snow within the layer k !
          !     [ kg/kg/s ]         by cumulus expels of condensates to the environment.!         
          ! --------------------------------------------------------------------------- !

           qrten(k) = frc_rasn * dwten(k)
           qsten(k) = frc_rasn * diten(k) 


          ! ------------------------------------------------------------------------ !
          ! 'slten(k)','qtten(k)'                                                    !
          !  Note that 'slflx(k)' and 'qtflx(k)' we have calculated already included !
          !  all the contributions of (1) expels of condensate (dwten(k), diten(k)), !
          !  (2) mass detrainment ( delta * umf * ( qtu - qt ) ), & (3) compensating !
          !  subsidence ( M * dqt / dz ). Thus 'slflx(k)' and 'qtflx(k)' we computed ! 
          !  is a hybrid turbulent flux containing one part of 'source' term - expel !
          !  of condensate. In order to calculate 'slten' and 'qtten', we should add !
          !  additional 'source' term, if any. If the expelled condensate falls down !
          !  across the base of the layer, it will be another sink (negative source) !
          !  term.  Note also that we included frictional heating terms in the below !
          !  calculation of 'slten'.                                                 !
          ! ------------------------------------------------------------------------ !
                   
           slten(k) = ( slflx(km1) - slflx(k) ) * g / dp0(k)
           if( k .eq. 1 ) then
              slten(k) = slten(k) - g / 4. / dp0(k) * (                            &
                                    uflx(k)*(uf(k+1) - uf(k) + u0(k+1) - u0(k)) +     & 
                                    vflx(k)*(vf(k+1) - vf(k) + v0(k+1) - v0(k)))
           elseif( k .ge. 2 .and. k .le. kpen-1 ) then
              slten(k) = slten(k) - g / 4. / dp0(k) * (                            &
                                    uflx(k)*(uf(k+1) - uf(k) + u0(k+1) - u0(k)) +     &
                                    uflx(k-1)*(uf(k) - uf(k-1) + u0(k) - u0(k-1)) +   &
                                    vflx(k)*(vf(k+1) - vf(k) + v0(k+1) - v0(k)) +     &
                                    vflx(k-1)*(vf(k) - vf(k-1) + v0(k) - v0(k-1)))
           elseif( k .eq. kpen ) then
              slten(k) = slten(k) - g / 4. / dp0(k) * (                            &
                                    uflx(k-1)*(uf(k) - uf(k-1) + u0(k) - u0(k-1)) +   &
                                    vflx(k-1)*(vf(k) - vf(k-1) + v0(k) - v0(k-1)))
           endif
           qtten(k) = ( qtflx(km1) - qtflx(k) ) * g / dp0(k)

          ! ---------------------------------------------------------------------------- !
          ! Compute condensate tendency, including reserved condensate                   !
          ! We assume that eventual detachment and detrainment occurs in kbup layer  due !
          ! to downdraft buoyancy sorting. In the layer above the kbup, only penetrative !
          ! entrainment exists. Penetrative entrained air is assumed not to contain any  !
          ! condensate.                                                                  !
          ! ---------------------------------------------------------------------------- !
  
          ! Compute in-cumulus condensate at the layer mid-point.

           if( k .lt. krel .or. k .gt. kpen ) then
              qlu_mid = 0.
              qiu_mid = 0.
              qlj     = 0.
              qij     = 0.
           elseif( k .eq. krel ) then 
              call conden(prel,thlu(krel-1),qtu(krel-1),thj,qvj,qlj,qij,qse,id_check)
              if( id_check .eq. 1 ) then
                  exit_conden(i) = 1.
                  id_exit = .true.
                  if (scverbose) then
                    call write_parallel('------- UW ShCu: exit, conden')
                  end if
                  go to 333
              endif
              qlubelow = qlj       
              qiubelow = qij       
              call conden(pifc0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check)
              if( id_check .eq. 1 ) then
                  exit_conden(i) = 1.
                  id_exit = .true.
                  if (scverbose) then
                    call write_parallel('------- UW ShCu: exit, conden')
                  end if
                  go to 333
              end if
              qlu_mid = 0.5 * ( qlubelow + qlj ) * ( prel - pifc0(k) )/( pifc0(k-1) - pifc0(k) )
              qiu_mid = 0.5 * ( qiubelow + qij ) * ( prel - pifc0(k) )/( pifc0(k-1) - pifc0(k) )
           elseif( k .eq. kpen ) then 
              call conden(pifc0(k-1)+ppen,thlu_top,qtu_top,thj,qvj,qlj,qij,qse,id_check)
              if( id_check .eq. 1 ) then
                  exit_conden(i) = 1.
                  id_exit = .true.
                  if (scverbose) then
                    call write_parallel('------- UW ShCu: exit, conden')
                  end if
                  go to 333
              end if
              qlu_mid = 0.5 * ( qlubelow + qlj ) * ( -ppen )        /( pifc0(k-1) - pifc0(k) )
              qiu_mid = 0.5 * ( qiubelow + qij ) * ( -ppen )        /( pifc0(k-1) - pifc0(k) )
              qlu_top = qlj
              qiu_top = qij
           else
              call conden(pifc0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check)
              if( id_check .eq. 1 ) then
                  exit_conden(i) = 1.
                  id_exit = .true.
                  if (scverbose) then
                    call write_parallel('------- UW ShCu: exit, conden')
                  end if
                  go to 333
              end if
              qlu_mid = 0.5 * ( qlubelow + qlj )
              qiu_mid = 0.5 * ( qiubelow + qij )
           endif
           qlubelow = qlj       
           qiubelow = qij

          ! 1. Non-precipitating portion of expelled condensate

           qc_l(k) = ( 1. - frc_rasn ) * dwten(k) ! [ kg/kg/s ]
           qc_i(k) = ( 1. - frc_rasn ) * diten(k) ! [ kg/kg/s ]

          ! 2. Detrained Condensate

           if( k .le. kbup ) then 
              qc_l(k) = qc_l(k) + g * 0.5 * ( umf(k-1) + umf(k) ) * fdr(k) * qlu_mid ! [ kg/kg/s ]
              qc_i(k) = qc_i(k) + g * 0.5 * ( umf(k-1) + umf(k) ) * fdr(k) * qiu_mid ! [ kg/kg/s ]
              qc_lm   =         - g * 0.5 * ( umf(k-1) + umf(k) ) * fdr(k) * ql0(k)  
              qc_im   =         - g * 0.5 * ( umf(k-1) + umf(k) ) * fdr(k) * qi0(k)
            ! Below 'nc_lm', 'nc_im' should be used only when frc_rasn = 1.
!              nc_lm   =         - g * 0.5 * ( umf(k-1) + umf(k) ) * fdr(k) * tr0(k,ixnumliq)  
!              nc_im   =         - g * 0.5 * ( umf(k-1) + umf(k) ) * fdr(k) * tr0(k,ixnumice)
           else
              qc_lm   = 0.
              qc_im   = 0.
              nc_lm   = 0.
              nc_im   = 0.
           endif

          ! 3. Detached Updraft 

           if( k .eq. kbup ) then
              qc_l(k) = qc_l(k) + g * umf(k) * qlj     / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
              qc_i(k) = qc_i(k) + g * umf(k) * qij     / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
              qc_lm   = qc_lm   - g * umf(k) * ql0(k)  / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
              qc_im   = qc_im   - g * umf(k) * qi0(k)  / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
!              nc_lm   = nc_lm   - g * umf(k) * tr0(k,ixnumliq)  / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
!              nc_im   = nc_im   - g * umf(k) * tr0(k,ixnumice)  / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
           endif 

          ! 4. Cumulative Penetrative entrainment detrained in the 'kbup' layer
          !    Explicitly compute the properties detrained penetrative entrained airs in k = kbup layer.

           if( k .eq. kbup ) then
              call conden(pmid0(k),thlu_emf(k),qtu_emf(k),thj,qvj,ql_emf_kbup,qi_emf_kbup,qse,id_check)
              if( id_check .eq. 1 ) then
                  id_exit = .true.
                  if (scverbose) then
                    call write_parallel('------- UW ShCu: exit, conden')
                  end if
                  go to 333
              endif
              if( ql_emf_kbup .gt. 0. ) then
!                  nl_emf_kbup = tru_emf(k,ixnumliq)
              else
                  nl_emf_kbup = 0.
              endif
              if( qi_emf_kbup .gt. 0. ) then
!                  ni_emf_kbup = tru_emf(k,ixnumice)
              else
                  ni_emf_kbup = 0.
              endif
              qc_lm   = qc_lm   - g * emf(k) * ( ql_emf_kbup - ql0(k) ) / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
              qc_im   = qc_im   - g * emf(k) * ( qi_emf_kbup - qi0(k) ) / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
!              nc_lm   = nc_lm   - g * emf(k) * ( nl_emf_kbup - tr0(k,ixnumliq) ) / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
!              nc_im   = nc_im   - g * emf(k) * ( ni_emf_kbup - tr0(k,ixnumice) ) / ( pifc0(k-1) - pifc0(k) ) ! [ kg/kg/s ]
           endif 

           qlten_det(k)   = qc_l(k) + qc_lm
           qiten_det(k)   = qc_i(k) + qc_im

          ! --------------------------------------------------------------------------------- !
          ! 'qlten(k)','qiten(k)','qvten(k)','sten(k)'                                        !
          !                                                                                   !
          ! --------------------------------------------------------------------------------- ! 

!           if( use_expconten ) then
!             if( use_unicondet ) then
!                  qc_l(k) = 0.
!                  qc_i(k) = 0. 
!                  qlten(k) = frc_rasn * dwten(k) + qlten_sink(k) + qlten_det
!                  qiten(k) = frc_rasn * diten(k) + qiten_sink(k) + qiten_det
!             else 
!                  qlten(k) = qc_l(k) + frc_rasn * dwten(k) + ( max( 0., ql0(k) + ( qc_lm + qlten_sink(k) ) * dt ) - ql0(k) ) / dt
!                  qiten(k) = qc_i(k) + frc_rasn * diten(k) + ( max( 0., qi0(k) + ( qc_im + qiten_sink(k) ) * dt ) - qi0(k) ) / dt
!                  trten(k,ixnumliq) = max( nc_lm + nlten_sink(k), - tr0(k,ixnumliq) / dt )
!                  trten(k,ixnumice) = max( nc_im + niten_sink(k), - tr0(k,ixnumice) / dt )
!             endif
!           else
!              if( use_unicondet ) then
!                  qc_l(k) = 0.
!                  qc_i(k) = 0. 
!              endif                      
!              qlten(k) = dwten(k) + ( qtten(k) - dwten(k) - diten(k) ) * ( ql0(k) / qt0(k) )
!              qiten(k) = diten(k) + ( qtten(k) - dwten(k) - diten(k) ) * ( qi0(k) / qt0(k) )
!           endif

           ! ----------------------------------------------------------------- !
           ! If combined sinks would result in negative liquid/ice, rescale    ! 
           ! such that liquid/ice is zero.                                     !
	   ! 
           ! ----------------------------------------------------------------- !

           if ( ((qc_lm+qlten_sink(k))*dt+ql0(k)).lt.0. ) then
              totsink = qc_lm+qlten_sink(k)
              if (totsink.ne.0.) then
                qc_lm = -(ql0(k)/dt) * qc_lm/totsink
                qlten_sink(k) = -(ql0(k)/dt) * qlten_sink(k)/totsink
                qlten_det(k) = qc_l(k) + qc_lm
              end if
           end if
           if ( ((qc_im+qiten_sink(k))*dt+qi0(k)).lt.0. ) then
              totsink = qc_im+qiten_sink(k)
              if (totsink.ne.0.) then
                qc_im = -(qi0(k)/dt) * qc_im/totsink
                qiten_sink(k) = -(qi0(k)/dt) * qiten_sink(k)/totsink
                qiten_det(k) = qc_i(k) + qc_im
              end if
           end if

           qlten(k) = qrten(k) + qlten_sink(k) + qlten_det(k)
           qiten(k) = qsten(k) + qiten_sink(k) + qiten_det(k)

           qvten(k) = qtten(k) - qlten(k) - qiten(k)
           sten(k)  = slten(k) + xlv * qlten(k) + xls * qiten(k)

!           qldet(k) = qc_l(k)
!           qidet(k) = qc_i(k)
           qc(k)  =  qc_l(k) +  qc_i(k)   


           qlten(k) = qlten(k) - qrten(k)
           qiten(k) = qiten(k) - qsten(k)
           qtten(k) = qlten(k) + qiten(k) + qvten(k)
!           if( ( qv0(k) + qvten(k)*dt ) .lt. 0.0 .or. &
!              ( ql0(k) + qlten(k)*dt ) .lt. 0.0 .or. &
!              ( qi0(k) + qiten(k)*dt ) .lt. 0.0 ) then
!               limit_negcon(i) = 1.
!           end if
           slten(k) = sten(k) - xlv*qlten(k) - xls*qiten(k)
           slten(k) = slten(k) + xlv * qrten(k) + xls * qsten(k)         
           sten(k)  = slten(k) + xlv * qlten(k) + xls * qiten(k)

         end do


       ! ---------------------------------------------------------------- !
       ! Now treats the 'evaporation' and 'melting' of rain ( qrten ) and ! 
       ! snow ( qsten ) during falling process.                           !
       ! ---------------------------------------------------------------- !

!         do k = k0, 1, -1  ! 'k' is a layer index : 'k0'('1') is the top ('bottom') layer
          
          ! ---------------------------------- !
          ! Calculate thermodynamic tendencies !
          ! --------------------------------------------------------------------------- !
          ! Note that equivalently, we can write tendency formula of 'sten' and 'slten' !
          ! by 'sten(k)  = sten(k) - xlv*evprain  - xls*evpsnow - (xls-xlv)*snowmlt' &  !
          !    'slten(k) = sten(k) - xlv*qlten(k) - xls*qiten(k)'.                      !
          ! The above formula is equivalent to the below formula. However below formula !
          ! is preferred since we have already imposed explicit constraint on 'ntraprd' !
          ! and 'ntsnprd' in case that flxrain(k-1) < 0 & flxsnow(k-1) < 0.          !
          ! Note : In future, I can elborate the limiting of 'qlten','qvten','qiten'    !
          !        such that that energy and moisture conservation error is completely  !
          !        suppressed.                                                          !
          ! Re-storation to the positive condensate will be performed later below       !
          ! --------------------------------------------------------------------------- !

!         end do

       ! --------------------------------------------------------------- !
       ! Prevent the onset-of negative condensate at the next time step  !
       ! Potentially, this block can be moved just in front of the above !
       ! block.                                                          ! 
       ! --------------------------------------------------------------- !

       ! Modification : I should check whether this 'positive_moisture_single' routine is
       !                consistent with the one used in UW PBL and cloud macrophysics schemes.
       ! Modification : Below may overestimate resulting 'ql, qi' if we use the new 'qc_l', 'qc_i'
       !                in combination with the original computation of qlten, qiten. However,
       !                if we use new 'qlten,qiten', there is no problem.

         qv0_star(:k0) = qv0(:k0) + qvten(:k0) * dt
         ql0_star(:k0) = ql0(:k0) + qlten(:k0) * dt
         qi0_star(:k0) = qi0(:k0) + qiten(:k0) * dt
         s0_star(:k0)  =  s0(:k0) +  sten(:k0) * dt
         call positive_moisture_single( xlv, xls, k0, dt, qmin(1), qmin(ixcldliq), qmin(ixcldice), &
              dp0, qv0_star, ql0_star, qi0_star, s0_star, qvten, qlten, qiten, sten )
         qtten(:k0)    = qvten(:k0) + qlten(:k0) + qiten(:k0)
         slten(:k0)    = sten(:k0)  - xlv * qlten(:k0) - xls * qiten(:k0)

       ! --------------------- !
       ! Tendencies of tracers !
       ! --------------------- !

       if (dotransport.eq.1) then
       do m = 1, ncnst

!         if( m .ne. ixnumliq .and. m .ne. ixnumice ) then

           trmin = 0. !qmin(m)
!#ifdef MODAL_AERO
!           do mm = 1, ntot_amode
!             if( m .eq. numptr_amode(mm) ) then
!                 trmin = 1.e-5
!                 goto 55
!             endif              
!           enddo
!        55 continue
!#endif 
           trflx_d(0:k0) = 0.
           trflx_u(0:k0) = 0.           
           do k = 1, k0-1
!             if( cnst_get_type_byind(m) .eq. 'wet' ) then
                 pdelx = dp0(k)
!             else
!                 pdelx = dpdry0(k)
!             endif
             km1 = k - 1
             dum = ( tr0(k,m) - trmin ) *  pdelx / g / dt + trflx(km1,m) - trflx(k,m) + trflx_d(km1)
             trflx_d(k) = min( 0., dum )
           enddo
           do k = k0, 2, -1
!             if( cnst_get_type_byind(m) .eq. 'wet' ) then
                 pdelx = dp0(k)
!             else
!                 pdelx = dpdry0(k)
!             endif
             km1 = k - 1
             dum = ( tr0(k,m) - trmin ) * pdelx / g / dt + trflx(km1,m) - trflx(k,m) + &
                                                           trflx_d(km1) - trflx_d(k) - trflx_u(k) 
             trflx_u(km1) = max( 0., -dum ) 
           enddo
           do k = 1, k0
!             if( cnst_get_type_byind(m) .eq. 'wet' ) then
                 pdelx = dp0(k)
!             else
!                 pdelx = dpdry0(k)
!             endif
             km1 = k - 1
           ! Check : I should re-check whether '_u', '_d' are correctly ordered in 
           !         the below tendency computation.
             trten(k,m) = ( trflx(km1,m) - trflx(k,m) + & 
                            trflx_d(km1) - trflx_d(k) + &
                            trflx_u(km1) - trflx_u(k) ) * g / pdelx
           enddo

!         endif

        enddo
        endif

       ! ---------------------------------------------------------------- !
       ! Cumpute default diagnostic outputs                               !
       ! Note that since 'qtu(krel-1:kpen-1)' & 'thlu(krel-1:kpen-1)' has !
       ! been adjusted after detraining cloud condensate into environment ! 
       ! during cumulus updraft motion,  below calculations will  exactly !
       ! reproduce in-cloud properties as shown in the output analysis.   !
       ! ---------------------------------------------------------------- ! 
 
       call conden(prel,thlu(krel-1),qtu(krel-1),thj,qvj,qlj,qij,qse,id_check)
       if( id_check .eq. 1 ) then
           exit_conden(i) = 1.
           id_exit = .true.
           if (scverbose) then
             call write_parallel('------- UW ShCu: exit, conden')
           end if
           go to 333
       end if
       qcubelow = qlj + qij
       qlubelow = qlj       
       qiubelow = qij       
       rcwp     = 0.
       rlwp     = 0.
       riwp     = 0.

       ! --------------------------------------------------------------------- !
       ! In the below calculations, I explicitly considered cloud base ( LCL ) !
       ! and cloud top height ( pifc0(kpen-1) + ppen )                           !
       ! ----------------------------------------------------------------------! 
       do k = krel, kpen ! This is a layer index
          ! ------------------------------------------------------------------ ! 
          ! Calculate cumulus condensate at the upper interface of each layer. !
          ! Note 'ppen < 0' and at 'k=kpen' layer, I used 'thlu_top'&'qtu_top' !
          ! which explicitly considered zero or non-zero 'fer(kpen)'.          !
          ! ------------------------------------------------------------------ ! 
          if( k .eq. kpen ) then 
              call conden(pifc0(k-1)+ppen,thlu_top,qtu_top,thj,qvj,qlj,qij,qse,id_check)
          else
              call conden(pifc0(k),thlu(k),qtu(k),thj,qvj,qlj,qij,qse,id_check)
          endif
          if( id_check .eq. 1 ) then
              exit_conden(i) = 1.
              id_exit = .true.
              if (scverbose) then
                call write_parallel('------- UW ShCu: exit, conden')
              end if
              go to 333
          end if
          ! ---------------------------------------------------------------- !
          ! Calculate in-cloud mean LWC ( qlu(k) ), IWC ( qiu(k) ),  & layer !
          ! mean cumulus fraction ( cufrc(k) ),  vertically-integrated layer !
          ! mean LWP and IWP. Expel some of in-cloud condensate at the upper !
          ! interface if it is largr than criqc. Note cumulus cloud fraction !
          ! is assumed to be twice of core updraft fractional area. Thus LWP !
          ! and IWP will be twice of actual value coming from our scheme.    !
          ! ---------------------------------------------------------------- !
          qcu(k)   = 0.5 * ( qcubelow + qlj + qij )
          qlu(k)   = 0.5 * ( qlubelow + qlj )
          qiu(k)   = 0.5 * ( qiubelow + qij )
          cufrc(k) = ( ufrc(k-1) + ufrc(k) )
          if( k .eq. krel ) then
              cufrc(k) = ( ufrclcl + ufrc(k) )*( prel - pifc0(k) )/( pifc0(k-1) - pifc0(k) )
          else if( k .eq. kpen ) then
              cufrc(k) = ( ufrc(k-1) + 0. )*( -ppen )        /( pifc0(k-1) - pifc0(k) )
              if( (qlj + qij) .gt. criqc ) then           
                   qcu(k) = 0.5 * ( qcubelow + criqc )
                   qlu(k) = 0.5 * ( qlubelow + criqc * qlj / ( qlj + qij ) )
                   qiu(k) = 0.5 * ( qiubelow + criqc * qij / ( qlj + qij ) )
              endif
          endif  
          rcwp = rcwp + ( qlu(k) + qiu(k) ) * ( pifc0(k-1) - pifc0(k) ) / g * cufrc(k)
          rlwp = rlwp +   qlu(k)            * ( pifc0(k-1) - pifc0(k) ) / g * cufrc(k)
          riwp = riwp +   qiu(k)            * ( pifc0(k-1) - pifc0(k) ) / g * cufrc(k)
          qcubelow = qlj + qij
          qlubelow = qlj
          qiubelow = qij
       end do
       ! ------------------------------------ !      
       ! Cloud top and base interface indices !
       ! ------------------------------------ !
       cnt = real( kpen )
       cnb = real( krel - 1 )

       ! ------------------------------------------------------------------------- !
       ! End of formal calculation. Below blocks are for implicit CIN calculations ! 
       ! with re-initialization and save variables at iter_cin = 1.             !
       ! ------------------------------------------------------------------------- !
       
       ! --------------------------------------------------------------- !
       ! Adjust the original input profiles for implicit CIN calculation !
       ! --------------------------------------------------------------- !

       if( iter .ne. iter_cin ) then 

          ! ------------------------------------------------------------------- !
          ! Save the output from "iter_cin = 1"                                 !
          ! These output will be writed-out if "iter_cin = 1" was not performed !
          ! for some reasons.                                                   !
          ! ------------------------------------------------------------------- !

           qv0_s(:k0)           = qv0(:k0) + qvten(:k0) * dt
           ql0_s(:k0)           = ql0(:k0) + qlten(:k0) * dt
           qi0_s(:k0)           = qi0(:k0) + qiten(:k0) * dt
           s0_s(:k0)            = s0(:k0)  +  sten(:k0) * dt 
           u0_s(:k0)            = u0(:k0)  +  uten(:k0) * dt
           v0_s(:k0)            = v0(:k0)  +  vten(:k0) * dt 
!           qt0_s(:k0)           = qv0_s(:k0) + ql0_s(:k0) + qi0_s(:k0)
           t0_s(:k0)            = t0(:k0)  +  sten(:k0) * dt / cp
           if (dotransport.eq.1) then
           do m = 1, ncnst
             tr0_s(:k0,m)       = tr0(:k0,m) + trten(:k0,m) * dt
           enddo
           endif

           umf_s(0:k0)          = umf(0:k0)
           dcm_s(:k0)           = dcm(:k0)
           qvten_s(:k0)         = qvten(:k0)
           qlten_s(:k0)         = qlten(:k0)  
           qiten_s(:k0)         = qiten(:k0)
           sten_s(:k0)          = sten(:k0)
           uten_s(:k0)          = uten(:k0)  
           vten_s(:k0)          = vten(:k0)
           qrten_s(:k0)         = qrten(:k0)
           qsten_s(:k0)         = qsten(:k0)  
           cush_s               = cush
           cufrc_s(:k0)         = cufrc(:k0)  
           slflx_s(0:k0)        = slflx(0:k0)  
           qtflx_s(0:k0)        = qtflx(0:k0)  
           uflx_s(0:k0)         = uflx(0:k0)  
           vflx_s(0:k0)         = vflx(0:k0)  
           qcu_s(:k0)           = qcu(:k0)  
           qlu_s(:k0)           = qlu(:k0)  
           qiu_s(:k0)           = qiu(:k0)  
           fer_s(:k0)           = fer(:k0)  
           fdr_s(:k0)           = fdr(:k0)  
           xc_s(:k0)            = xco(:k0)
           cin_s                 = cin
           cinlcl_s              = cinlcl
           cbmf_s                = cbmf
           qc_s(:k0)            = qc(:k0)
           qldet_s(:k0)         = qlten_det(:k0)
           qidet_s(:k0)         = qiten_det(:k0)
           qlsub_s(:k0)         = qlten_sink(:k0)
           qisub_s(:k0)         = qiten_sink(:k0)

!           slten_s(:k0)         = slten(:k0)
           ufrc_s(0:k0)         = ufrc(0:k0) 


#ifdef UWDIAG         
           cnt_s                = cnt
           cnb_s                = cnb
           qtten_s(:k0)         = qtten(:k0)
           ufrcinvbase_s        = ufrcinvbase
           ufrclcl_s            = ufrclcl 
           winvbase_s           = winvbase
           wlcl_s               = wlcl
           plcl_s               = plcl
           pinv_s               = pifc0(kinv-1)
           plfc_s               = plfc        
           prel_s               = prel        
           pbup_s               = pifc0(kbup)
           ppen_s               = pifc0(kpen-1) + ppen        
           qtsrc_s              = qtsrc
           thlsrc_s             = thlsrc
           thvlsrc_s            = thvlsrc
           emfkbup_s            = emf(kbup)
           cbmflimit_s          = cbmflimit
           tkeavg_s             = tkeavg
           zinv_s               = zifc0(kinv-1)
           rcwp_s               = rcwp
           rlwp_s               = rlwp
           riwp_s               = riwp
           wu_s(0:k0)           = wu(0:k0)
           qtu_s(0:k0)          = qtu(0:k0)
           thlu_s(0:k0)         = thlu(0:k0)
           thvu_s(0:k0)         = thvu(0:k0)
           uu_s(0:k0)           = uu(0:k0)
           vu_s(0:k0)           = vu(0:k0)
           qtu_emf_s(0:k0)      = qtu_emf(0:k0)
           thlu_emf_s(0:k0)     = thlu_emf(0:k0)
           uu_emf_s(0:k0)       = uu_emf(0:k0)
           vu_emf_s(0:k0)       = vu_emf(0:k0)
           uemf_s(0:k0)         = uemf(0:k0)

           dwten_s(:k0)         = dwten(:k0)
           diten_s(:k0)         = diten(:k0)

           excessu_arr_s(:k0)   = excessu_arr(:k0)
           excess0_arr_s(:k0)   = excess0_arr(:k0)
           xc_arr_s(:k0)        = xc_arr(:k0)
           aquad_arr_s(:k0)     = aquad_arr(:k0)
           bquad_arr_s(:k0)     = bquad_arr(:k0)
           cquad_arr_s(:k0)     = cquad_arr(:k0)
           bogbot_arr_s(:k0)    = bogbot_arr(:k0)
           bogtop_arr_s(:k0)    = bogtop_arr(:k0)

           if (dotransport.eq.1) then
           do m = 1, ncnst
             trten_s(:k0,m)    = trten(:k0,m)
             trflx_s(0:k0,m)   = trflx(0:k0,m)
             tru_s(0:k0,m)     = tru(0:k0,m)
             tru_emf_s(0:k0,m) = tru_emf(0:k0,m)
           enddo
           endif
#endif



          ! ----------------------------------------------------------------------------- ! 
          ! Recalculate environmental variables for new cin calculation at "iter_cin = 2" ! 
          ! using the updated state variables. Perform only for variables necessary  for  !
          ! the new cin calculation.                                                      !
          ! ----------------------------------------------------------------------------- !
          
           qv0(:k0)   = qv0_s(:k0)
           ql0(:k0)   = ql0_s(:k0)
           qi0(:k0)   = qi0_s(:k0)
           s0(:k0)    = s0_s(:k0)
           t0(:k0)    = t0_s(:k0)
      
           qt0(:k0)   = (qv0(:k0) + ql0(:k0) + qi0(:k0))
           thl0(:k0)  = (t0(:k0) - xlv*ql0(:k0)/cp - xls*qi0(:k0)/cp)/exnmid0(:k0)
           thvl0(:k0) = (1. + zvir*qt0(:k0))*thl0(:k0)

           ssthl0      = slope(k0,thl0,pmid0) ! Dimension of ssthl0(:k0) is implicit
           ssqt0       = slope(k0,qt0 ,pmid0)
           ssu0        = slope(k0,u0  ,pmid0)
           ssv0        = slope(k0,v0  ,pmid0)
           if (dotransport.eq.1) then
           do m = 1, ncnst
             sstr0(:k0,m) = slope(k0,tr0(:k0,m),pmid0)
           enddo
           endif

           do k = 1, k0

             thl0bot = thl0(k) + ssthl0(k) * ( pifc0(k-1) - pmid0(k) )
             qt0bot  = qt0(k)  + ssqt0(k)  * ( pifc0(k-1) - pmid0(k) )
             call conden(pifc0(k-1),thl0bot,qt0bot,thj,qvj,qlj,qij,qse,id_check)
             if( id_check .eq. 1 ) then
                 exit_conden(i) = 1.
                 id_exit = .true.
                 if (scverbose) then
                   call write_parallel('------- UW ShCu: exit, conden')
                 end if
                 go to 333
             end if
             thv0bot(k)  = thj * ( 1. + zvir*qvj - qlj - qij )
             thvl0bot(k) = thl0bot * ( 1. + zvir*qt0bot )
          
             thl0top = thl0(k) + ssthl0(k) * ( pifc0(k) - pmid0(k) )
             qt0top  =  qt0(k) + ssqt0(k)  * ( pifc0(k) - pmid0(k) )
             call conden(pifc0(k),thl0top,qt0top,thj,qvj,qlj,qij,qse,id_check)
             if( id_check .eq. 1 ) then
                 exit_conden(i) = 1.
                 id_exit = .true.
                 if (scverbose) then
                   call write_parallel('------- UW ShCu: exit, conden')
                 end if
                 go to 333
             end if
             thv0top(k)  = thj * ( 1. + zvir*qvj - qlj - qij )
             thvl0top(k) = thl0top * ( 1. + zvir*qt0top )

           end do

         endif               ! End of 'if(iter .ne. iter_cin)' if sentence. 

       end do  ! iter loop



     ! ----------------------- !
     ! Update Output Variables !
     ! ----------------------- !

     umf_out(i,0:k0)             = umf(0:k0)
     umf_out(i,0:kinv-2)         = uemf(0:kinv-2)

!     umf_out(i,0:kinv-2)         = uemf(0:kinv-2)
     dcm_out(i,:k0)              = dcm(:k0)
!the indices are not reversed, these variables go into compute_mcshallow_inv
     qvten_out(i,:k0)            = qvten(:k0)
     qlten_out(i,:k0)            = qlten(:k0)
     qiten_out(i,:k0)            = qiten(:k0)
     sten_out(i,:k0)             = sten(:k0)
     uten_out(i,:k0)             = uten(:k0)
     vten_out(i,:k0)             = vten(:k0)
     qrten_out(i,:k0)            = qrten(:k0)
     qsten_out(i,:k0)            = qsten(:k0)
     cufrc_out(i,:k0)            = cufrc(:k0)
     cush_inout(i)               = cush
     qldet_out(i,:k0)            = qlten_det(:k0)
     qidet_out(i,:k0)            = qiten_det(:k0)
     qlsub_out(i,:k0)            = qlten_sink(:k0)
     qisub_out(i,:k0)            = qiten_sink(:k0)
     ndrop_out(i,:k0)            = qlten_det(:k0)/(4188.787*rdrop**3)
!     ndrop_out(i,:k0)            = qlten_det(:k0)/(4.19e-12) !(1.15e-11) ! /drop mass
     nice_out(i,:k0)             = qiten_det(:k0)/(3.0e-10) ! /crystal mass
     qtflx_out(i,0:k0)           = qtflx(0:k0)
     slflx_out(i,0:k0)           = slflx(0:k0)
     uflx_out(i,0:k0)            = uflx(0:k0)
     vflx_out(i,0:k0)            = vflx(0:k0)

     if (dotransport.eq.1) then
     do m = 1, ncnst
        tr0_inout(i,:k0,m)      = tr0_inout(i,:k0,m) + trten(:k0,m) * dt
     enddo
     endif
  
     ! ------------------------------------------------- !
     ! Below are specific diagnostic output for detailed !
     ! analysis of cumulus scheme                        !
     ! ------------------------------------------------- !

     fer_out(i,1:kpen)          = fer(:kpen)  
     fdr_out(i,1:kpen)          = fdr(:kpen)  

#ifdef UWDIAG
     cldhgt_out(i)               = cldhgt
     cbmf_out(i)                 = cbmf
     cnt_out(i)                  = cnt
     cnb_out(i)                  = cnb
     qcu_out(i,:k0)              = qcu(:k0)
     qlu_out(i,:k0)              = qlu(:k0)
     qiu_out(i,:k0)              = qiu(:k0)
     qc_out(i,:k0)               = qc(:k0)
     xc_out(i,1:k0)           = xco(:k0)
     cinh_out(i)              = cin
     cinlclh_out(i)           = cinlcl
!     qtten_out(i,1:k0)        = qtten(:k0)
!     slten_out(i,1:k0)        = slten(:k0)
!     ufrc_out(i,0:k0)         = ufrc(0:k0)
!     uflx_out(i,0:k0)         = uflx(0:k0)  
!     vflx_out(i,0:k0)         = vflx(0:k0)  
     
     ufrcinvbase_out(i)           = ufrcinvbase
     ufrclcl_out(i)               = ufrclcl 
     winvbase_out(i)              = winvbase
     wlcl_out(i)                  = wlcl
     plcl_out(i)                  = plcl
     pinv_out(i)                  = pifc0(kinv-1)
     plfc_out(i)                  = plfc    
     prel_out(i)                  = prel    
     pbup_out(i)                  = pifc0(kbup)        
     ppen_out(i)                  = pifc0(kpen-1) + ppen            
     qtsrc_out(i)                 = qtsrc
     thlsrc_out(i)                = thlsrc
     thvlsrc_out(i)               = thvlsrc
     emfkbup_out(i)               = emf(kbup)
     cbmflimit_out(i)             = cbmflimit
     tkeavg_out(i)                = tkeavg
     zinv_out(i)                  = zifc0(kinv-1)
     rcwp_out(i)                  = rcwp
     rlwp_out(i)                  = rlwp
     riwp_out(i)                  = riwp

     wu_out(i,0:k0)           = wu(0:k0)
     qtu_out(i,0:k0)          = qtu(0:k0)
     thlu_out(i,0:k0)         = thlu(0:k0)
     thvu_out(i,0:k0)         = thvu(0:k0)
     uu_out(i,0:k0)           = uu(0:k0)
     vu_out(i,0:k0)           = vu(0:k0)
     qtu_emf_out(i,0:k0)      = qtu_emf(0:k0)
     thlu_emf_out(i,0:k0)     = thlu_emf(0:k0)
     uu_emf_out(i,0:k0)       = uu_emf(0:k0)
     vu_emf_out(i,0:k0)       = vu_emf(0:k0)
     uemf_out(i,0:k0)         = uemf(0:k0)

     dwten_out(i,1:k0)        = dwten(:k0)
     diten_out(i,1:k0)        = diten(:k0)

     excessu_arr_out(i,1:k0)  = excessu_arr(:k0)
     excess0_arr_out(i,1:k0)  = excess0_arr(:k0)
     xc_arr_out(i,1:k0)       = xc_arr(:k0)
     aquad_arr_out(i,1:k0)    = aquad_arr(:k0)
     bquad_arr_out(i,1:k0)    = bquad_arr(:k0)
     cquad_arr_out(i,1:k0)    = cquad_arr(:k0)
     bogbot_arr_out(i,1:k0)   = bogbot_arr(:k0)
     bogtop_arr_out(i,1:k0)   = bogtop_arr(:k0)

      if (dotransport.eq.1) then
      do m = 1, ncnst
        trten_out(i,:k0,m)    = trten(:k0,m)
        trflx_out(i,0:k0,m)   = trflx(0:k0,m)  
        tru_out(i,0:k0,m)     = tru(0:k0,m)
        tru_emf_out(i,0:k0,m) = tru_emf(0:k0,m)
      enddo
      endif
#endif

333    if (id_exit) then

         exit_uwcu(i) = 1.
         if (scverbose) then
           call write_parallel('------- UW ShCu: Exited!')
         end if

     ! --------------------------------------------------------------------- !
     ! Initialize output variables when cumulus convection was not performed.!
     ! --------------------------------------------------------------------- !
     
     umf_out(i,0:k0)             = 0.   
     dcm_out(i,:k0)              = 0.   
     qvten_out(i,:k0)            = 0.
     qlten_out(i,:k0)            = 0.
     qiten_out(i,:k0)            = 0.
     sten_out(i,:k0)             = 0.
     uten_out(i,:k0)             = 0.
     vten_out(i,:k0)             = 0.
     qrten_out(i,:k0)            = 0.
     qsten_out(i,:k0)            = 0.
     cufrc_out(i,:k0)            = 0.
     cush_inout(i)               = -1.
     qldet_out(i,:k0)            = 0.
     qidet_out(i,:k0)            = 0.
     qtflx_out(i,0:k0)           = 0.
     slflx_out(i,0:k0)           = 0.
     uflx_out(i,0:k0)            = 0.
     vflx_out(i,0:k0)            = 0.

     fer_out(i,1:k0)             = MAPL_UNDEF
     fdr_out(i,1:k0)             = MAPL_UNDEF

#ifdef UWDIAG
     cbmf_out(i)                 = 0.   
     cnt_out(i)                  = 1.
     cnb_out(i)                  = real(k0)
     qcu_out(i,:k0)              = 0.
     qlu_out(i,:k0)              = 0.
     qiu_out(i,:k0)              = 0.
     qc_out(i,:k0)               = 0.
     xc_out(i,1:k0)              = MAPL_UNDEF
     cinh_out(i)                 = cin 
     cinlclh_out(i)              = cinlcl
!     qtten_out(i,k0:1:-1)        = 0.
!     slten_out(i,k0:1:-1)        = 0.
!     ufrc_out(i,k0:0:-1)         = 0.
!     uflx_out(i,k0:0:-1)         = 0.  
!     vflx_out(i,k0:0:-1)         = 0.  

     ufrcinvbase_out(i)           = 0. 
     ufrclcl_out(i)               = 0. 
     winvbase_out(i)              = 0.    
     wlcl_out(i)                  = MAPL_UNDEF    
     plcl_out(i)                  = MAPL_UNDEF
     pinv_out(i)                  = MAPL_UNDEF
     prel_out(i)                  = MAPL_UNDEF
     plfc_out(i)                  = MAPL_UNDEF
     pbup_out(i)                  = MAPL_UNDEF
     ppen_out(i)                  = MAPL_UNDEF
     qtsrc_out(i)                 = MAPL_UNDEF 
     thlsrc_out(i)                = MAPL_UNDEF    
     thvlsrc_out(i)               = MAPL_UNDEF
     emfkbup_out(i)               = 0.
     cbmflimit_out(i)             = 0.    
     tkeavg_out(i)                = tkeavg    
     zinv_out(i)                  = 0.    
     rcwp_out(i)                  = 0.    
     rlwp_out(i)                  = 0.    
     riwp_out(i)                  = 0.    

     wu_out(i,k0:0:-1)           = MAPL_UNDEF
     qtu_out(i,k0:0:-1)          = MAPL_UNDEF
     thlu_out(i,k0:0:-1)         = MAPL_UNDEF 
     thvu_out(i,k0:0:-1)         = MAPL_UNDEF 
     uu_out(i,k0:0:-1)           = MAPL_UNDEF
     vu_out(i,k0:0:-1)           = MAPL_UNDEF
     qtu_emf_out(i,k0:0:-1)      = MAPL_UNDEF
     thlu_emf_out(i,k0:0:-1)     = MAPL_UNDEF         
     uu_emf_out(i,k0:0:-1)       = MAPL_UNDEF  
     vu_emf_out(i,k0:0:-1)       = MAPL_UNDEF
     uemf_out(i,k0:0:-1)         = MAPL_UNDEF
   
     dwten_out(i,k0:1:-1)        = 0.    
     diten_out(i,k0:1:-1)        = 0.    

        excessu_arr_out(i,k0:1:-1)  = 0.    
        excess0_arr_out(i,k0:1:-1)  = 0.    
        xc_arr_out(i,k0:1:-1)       = 0.    
        aquad_arr_out(i,k0:1:-1)    = 0.    
        bquad_arr_out(i,k0:1:-1)    = 0.    
        cquad_arr_out(i,k0:1:-1)    = 0.    
        bogbot_arr_out(i,k0:1:-1)   = 0.    
        bogtop_arr_out(i,k0:1:-1)   = 0.    

        if (dotransport.eq.1) then
        do m = 1, ncnst
          trten_out(i,:k0,m)       = 0.
          trflx_out(i,k0:0:-1,m)   = 0.  
          tru_out(i,k0:0:-1,m)     = 0.
          tru_emf_out(i,k0:0:-1,m) = 0.
        enddo
        endif
#endif

       end if

     end do ! column i loop

     return

   end subroutine compute_uwshcu


   function slope(k0,field,p0)

      integer, intent(in):: k0
      real             :: slope(k0)
      real, intent(in) :: field(k0)
      real, intent(in) :: p0(k0)
      real             :: below
      real             :: above
      integer            :: k

      below = ( field(2) - field(1) ) / ( p0(2) - p0(1) )
      do k = 2, k0
         above = ( field(k) - field(k-1) ) / ( p0(k) - p0(k-1) )  
         if ( above .gt. 0. ) then
            slope(k-1) = max(0.,min(above,below))
         else
            slope(k-1) = min(0.,max(above,below))
         end if
         below = above
      end do
      slope(k0) = slope(k0-1)

      return
   end function slope


  function qsinvert(qt,thl,ps_in)
  ! ----------------------------------------------------------------- !
  ! Function calculating saturation pressure ps (or pLCL) from qt and !
  ! thl ( liquid potential temperature,  NOT liquid virtual potential ! 
  ! temperature) by inverting Bolton formula. I should check later if !
  ! current use of 'leff' instead of 'xlv' here is reasonable or not. !
  ! ----------------------------------------------------------------- !

     implicit none

     real          :: ps_in
     real          :: qsinvert    
     real,intent(in) :: qt, thl
     real*8          :: Pis, err, dlnqsdT, dTdPis, dqsdT
     real*8          :: dPisdps, dlnqsdps, derrdps, dps 
     real*8          :: rhi, TLCL, PiLCL, psmin, dpsmax
     real*8            :: Ti,qs,Ts,ps
     real            :: Tgeos,Qgeos,Pgeos
     integer         :: i
     real*8          :: es                     ! saturation vapor pressure
     real*8          :: gam                    ! (L/cp)*dqs/dT
     real*8          :: leff, nu

     psmin  = 10000._r8 ! Default saturation pressure [Pa] if iteration does not converge
     dpsmax = 1._r8           ! Tolerance [Pa] for convergence of iteration

     ! ------------------------------------ !
     ! Calculate best initial guess of pLCL !
     ! ------------------------------------ !

     Ti       =  thl*(ps_in/p00)**rovcp
     Tgeos    = Ti
     Pgeos    = ps_in
     qs       = dble( GEOS_QSAT(Tgeos,Pgeos/100.) )
     es       = ps_in * qs  / ( ep2 + (1._r8-ep2)*qs )
     rhi      = qt/qs
     if( rhi .le. 0.01_r8 ) then
!        if (scverbose) then
          call write_parallel('Source air is too dry and pLCL is set to psmin in uwshcu.F90')
!        end if
        qsinvert = psmin
        return
     end if

!     print *,'Ti,Rhi,thl=',Ti,rhi,thl
     TLCL     =  55._r8 + 1._r8/(1._r8/(Ti-55._r8)-log(rhi)/2840._r8) ! Bolton's formula. MWR.1980.Eq.(22)
     PiLCL    =  TLCL/thl
     ps       =  p00*(PiLCL)**(1._r8/rovcp)

     do i = 1, 10
       Pis      =  (ps/p00)**rovcp   ! Exner function
       Ts       =  thl*Pis
       Tgeos    = Ts
       Pgeos    = ps
       qs       = dble( GEOS_QSAT(Tgeos,Pgeos/100.) )
       Qgeos    = qs
       dqsdT    = dble( GEOS_DQSAT(Tgeos, Pgeos/100., QSAT=Qgeos ) )
       gam      = (xlv/cp)*dqsdT
       err      =  qt - qs
       nu       = 1.0 - fract_liq_f( real(Ts,4) )
!       nu       =  max(min((258._r8 - Ts)/20._r8,1.0_r8),0.0_r8)        
       leff     =  (1._r8 - nu)*xlv + nu*xls                  
       dlnqsdT  =  gam*(cp/leff)/qs
       dTdPis   =  thl
       dPisdps  =  rovcp*Pis/ps 
       dlnqsdps = -1._r8/(ps - (1. - ep2)*es)
       derrdps  = -qs*(dlnqsdT * dTdPis * dPisdps + dlnqsdps)
       if (derrdps.eq.0._r8) then
         print *,"QSINVERT: derrdps=0 !!!"
       end if
       dps      = -err/derrdps
       ps       =  ps + dps
       if( ps .lt. 0._r8 ) then
!           if (scverbose) then
!             call write_parallel('------- UW ShCu: pLCL iteration is negative and set to psmin')
!           end if
           qsinvert = psmin
           return    
       end if
       if( abs(dps) .le. dpsmax ) then
           qsinvert = ps
           return
       end if
     end do
!     if (scverbose) then
!       call write_parallel('------- UW ShCu: pLCL does not converge and is set to psmin')
!     end if
     qsinvert = psmin
     return
   end function qsinvert


   real function compute_alpha(del_CIN,ke)
   ! ------------------------------------------------ !
   ! Subroutine to compute proportionality factor for !
   ! implicit CIN calculation.                        !   
   ! ------------------------------------------------ !
     real   :: del_CIN, ke
     real*8 :: del_CIN8, ke8
     real*8 :: x0, x1

     integer  :: iteration

     x0 = 0._r8
     del_CIN8 = del_CIN
     ke8 = ke
     do iteration = 1, 10
        x1 = x0 - (exp(-x0*ke8*del_CIN8) - x0)/(-ke8*del_CIN8*exp(-x0*ke8*del_CIN8) - 1.)
        x0 = x1
     end do
     compute_alpha = x0

     return

   end function compute_alpha   


  real function compute_mumin2(mulcl,rmaxfrac,mulow)
  ! --------------------------------------------------------- !
  ! Subroutine to compute critical 'mu' (normalized CIN) such ! 
  ! that updraft fraction at the LCL is equal to 'rmaxfrac'.  !
  ! --------------------------------------------------------- !  
    real :: mulcl, rmaxfrac, mulow
    real*8 :: x0, x1, ex, ef, exf, f, fs
    integer  :: iteration

    x0 = mulow
    do iteration = 1, 10
       ex = exp(-x0**2)
       ef = erfc(x0)
       exf = ex/ef
       f  = 0.5_r8*exf**2 - 0.5_r8*(ex/2._r8/rmaxfrac)**2 - (mulcl*2.5066_r8/2._r8)**2
       fs = (2._r8*exf**2)*(exf/sqrt(MAPL_PI)-x0) + (0.5_r8*x0*ex**2)/(rmaxfrac**2)
       x1 = x0 - f/fs     
       x0 = x1
    end do
    compute_mumin2 = x0

    return

  end function compute_mumin2


  real function compute_ppen(wtwb,D,bogbot,bogtop,rho0j,dpen)
  ! ----------------------------------------------------------- !
  ! Subroutine to compute critical 'ppen[Pa]<0' ( pressure dis. !
  ! from 'pifc0(kpen-1)' to the cumulus top where cumulus updraft !
  ! vertical velocity is exactly zero ) by considering exact    !
  ! non-zero fer(kpen).                                         !  
  ! ----------------------------------------------------------- !  
    real :: wtwb, D, bogbot, bogtop, rho0j, dpen
    real*8 :: x0, x1, f, fs, SB, s00, aux
    integer  :: iteration

    ! Buoyancy slope
      SB = ( bogtop - bogbot ) / dpen
    ! Sign of slope, 'f' at x = 0
    ! If 's00>0', 'w' increases with height.
      s00 = bogbot / rho0j - D * wtwb

    if( D*dpen .lt. 1.e-4_r8 ) then
        if( s00 .ge. 0._r8 ) then
            x0 = dpen       
        else
            x0 = max(0._r8,min(dpen,-0.5_r8*wtwb/s00))
        endif
    else
        if( s00 .ge. 0._r8 ) then
            x0 = dpen
        else 
            x0 = 0._r8
        endif
        do iteration = 1, 5
           aux  = min(max(-2._r8*D*x0, -20.0), 20.0)
           
           f  = exp(aux)*(wtwb-(bogbot-SB/(2._r8*D))/(D*rho0j)) + &
                                 (SB*x0+bogbot-SB/(2._r8*D))/(D*rho0j)
           fs = -2._r8*D*exp(aux)*(wtwb-(bogbot-SB/(2._r8*D))/(D*rho0j)) + &
                                 (SB)/(D*rho0j)
!           if( fs .ge. 0._r8 ) then
!	     fs = max(fs, 1.e-8_r8)
!           else
!             fs = min(fs,-1.e-8_r8)
!           endif
           
           x1 = x0 - f/fs     
           x0 = x1
!           if (x0.gt.dpen*10._r8) then
!              print *,'compute_ppen: x0 limit reached'
!              x0 = 0._r8
!              exit
!           end if
        end do
    endif    

    compute_ppen = -max(0._r8,min(dpen,x0))

  end function compute_ppen


  subroutine fluxbelowinv(cbmf,ps0,mkx,kinv,dt,xsrc,xmean,xtopin,xbotin,xflx)   
  ! ------------------------------------------------------------------------- !
  ! Subroutine to calculate turbulent fluxes at and below 'kinv-1' interfaces.!
  ! Check in the main program such that input 'cbmf' should not be zero.      !  
  ! If the reconstructed inversion height does not go down below the 'kinv-1' !
  ! interface, then turbulent flux at 'kinv-1' interface  is simply a product !
  ! of 'cmbf' and 'qtsrc-xbot' where 'xbot' is the value at the top interface !
  ! of 'kinv-1' layer. This flux is linearly interpolated down to the surface !
  ! assuming turbulent fluxes at surface are zero. If reconstructed inversion !
  ! height goes down below the 'kinv-1' interface, subsidence warming &drying !
  ! measured by 'xtop-xbot', where  'xtop' is the value at the base interface !
  ! of 'kinv+1' layer, is added ONLY to the 'kinv-1' layer, using appropriate !
  ! mass weighting ( rpinv and rcbmf, or rr = rpinv / rcbmf ) between current !
  ! and next provisional time step. Also impose a limiter to enforce outliers !
  ! of thermodynamic variables in 'kinv' layer  to come back to normal values !
  ! at the next step.                                                         !
  ! ------------------------------------------------------------------------- !            
    integer,  intent(in)                     :: mkx, kinv 
    real, intent(in)                     :: cbmf, xsrc, xmean, xtopin, xbotin
    real, intent(in),  dimension(0:mkx)  :: ps0
    real, intent(out), dimension(0:mkx)  :: xflx  
    integer k
    real, intent(in) ::  dt
    real rcbmf, rpeff, dp, rr, pinv_eff, xtop, xbot, pinv, xtop_ori, xbot_ori

    xflx(0:mkx) = 0.
    dp = ps0(kinv-1) - ps0(kinv)    
    xbot = xbotin
    xtop = xtopin
   
    ! -------------------------------------- !
    ! Compute reconstructed inversion height !
    ! -------------------------------------- !
    xtop_ori = xtop
    xbot_ori = xbot
    rcbmf = ( cbmf * g * dt ) / dp                  ! Can be larger than 1 : 'OK'      

    if( xbot .ge. xtop ) then
        rpeff = ( xmean - xtop ) / max(  1.e-20, xbot - xtop ) 
    else
        rpeff = ( xmean - xtop ) / min( -1.e-20, xbot - xtop ) 
    endif 

    rpeff = min( max(0.,rpeff), 1. )          ! As of this, 0<= rpeff <= 1   
    if( rpeff .eq. 0. .or. rpeff .eq. 1. ) then
        xbot = xmean
        xtop = xmean
    endif
    ! Below two commented-out lines are the old code replacing the above 'if' block.   
    ! if(rpeff.eq.1) xbot = xmean
    ! if(rpeff.eq.0) xtop = xmean    
    rr       = rpeff / rcbmf
    pinv     = ps0(kinv-1) - rpeff * dp             ! "pinv" before detraining mass
    pinv_eff = ps0(kinv-1) + ( rcbmf - rpeff ) * dp ! Effective "pinv" after detraining mass
    ! ----------------------------------------------------------------------- !
    ! Compute turbulent fluxes.                                               !
    ! Below two cases exactly converges at 'kinv-1' interface when rr = 1. !
    ! ----------------------------------------------------------------------- !
    do k = 0, kinv - 1
       xflx(k) = cbmf * ( xsrc - xbot ) * ( ps0(0) - ps0(k) ) / ( ps0(0) - pinv )
    end do
    if( rr .le. 1. ) then
        xflx(kinv-1) =  xflx(kinv-1) - ( 1. - rr ) * cbmf * ( xtop_ori - xbot_ori )
    endif

    return
  end subroutine fluxbelowinv


  subroutine positive_moisture_single( xlv, xls, mkx, dt, qvmin, qlmin, qimin, dp, qv, ql, qi, s, qvten, qlten, qiten, sten )
  ! ------------------------------------------------------------------------------- !
  ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
  ! force them to be larger than minimum value by (1) condensating water vapor      !
  ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
  ! layer. '2._r8' is multiplied to the minimum values for safety.                  !
  ! Update final state variables and tendencies associated with this correction.    !
  ! If any condensation happens, update (s,t) too.                                  !
  ! Note that (qv,ql,qi,s) are final state variables after applying corresponding   !
  ! input tendencies and corrective tendencies                                      !
  ! ------------------------------------------------------------------------------- !
    implicit none
    integer,  intent(in)     :: mkx
    real, intent(in)     :: xlv, xls
    real, intent(in)     :: qvmin, qlmin, qimin
    real, intent(in)     :: dp(mkx)
    real, intent(in)       :: dt
    real, intent(inout)  :: qv(mkx), ql(mkx), qi(mkx), s(mkx)
    real, intent(inout)  :: qvten(mkx), qlten(mkx), qiten(mkx), sten(mkx)
    integer   k
    real*8 dql, dqi, dqv, sum, aa, dum 

    do k = mkx, 1, -1        ! From the top to the 1st (lowest) layer from the surface
       dql = max(0._r8,1._r8*qlmin-ql(k))
       dqi = max(0._r8,1._r8*qimin-qi(k))
       qlten(k) = qlten(k) +  dql/dt
       qiten(k) = qiten(k) +  dqi/dt
       qvten(k) = qvten(k) - (dql+dqi)/dt
       sten(k)  = sten(k)  + xlv * (dql/dt) + xls * (dqi/dt)
       ql(k)    = ql(k) +  dql
       qi(k)    = qi(k) +  dqi
       qv(k)    = qv(k) -  dql - dqi
       s(k)     = s(k)  +  xlv * dql + xls * dqi
       dqv      = max(0.,1.*qvmin-qv(k))
       qvten(k) = qvten(k) + dqv/dt
       qv(k)    = qv(k)   + dqv
       if( k .ne. 1 ) then 
           qv(k-1)    = qv(k-1)    - dqv*dp(k)/dp(k-1)
           qvten(k-1) = qvten(k-1) - dqv*dp(k)/dp(k-1)/dt
       endif
       qv(k) = max(qv(k),qvmin)
       ql(k) = max(ql(k),qlmin)
       qi(k) = max(qi(k),qimin)
    end do
    ! Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally 
    ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
    ! preserves column moisture. 
    if( dqv .gt. 1.e-20_r8 ) then
        sum = 0.
        do k = 1, mkx
           if( qv(k) .gt. 2._r8*qvmin ) sum = sum + qv(k)*dp(k)
        enddo
        aa = dqv*dp(1)/max(1.e-20_r8,sum)
        if( aa .lt. 0.5_r8 ) then
            do k = 1, mkx
               if( qv(k) .gt. 2._r8*qvmin ) then
                   dum      = aa*qv(k)
                   qv(k)    = qv(k) - dum
                   qvten(k) = qvten(k) - dum/dt
               endif
            enddo 
        else 
!           if (scverbose) then
!             call write_parallel( 'Full positive_moisture is impossible in uwshcu')
!           end if
        endif
    endif 

    return
  end subroutine positive_moisture_single


   subroutine conden(p,thl,qt,th,qv,ql,qi,rvls,id_check)

      implicit none

      real, intent(in)     :: p
      real, intent(in)     :: thl
      real, intent(in)     :: qt
      real, intent(out)    :: th
      real, intent(out)    :: qv
      real, intent(out)    :: ql
      real, intent(out)    :: qi
      real, intent(out)    :: rvls
      integer, intent(out) :: id_check

      ! Local
      real*8               :: tc
      real                 :: temps,ps
      real*8               :: leff,nu,qc
      integer              :: iteration
      real*8               :: qs    ! Saturation specific humidity

      tc = thl*exnerfn(p)

!      nu = max(min((258._r8-tc)/20._r8,1.0_r8),0.0_r8)    ! ice fraction of condensate
      nu = 1.0 - fract_liq_f( real(tc,4) )
      leff = (1._r8- nu)*xlv + nu*xls    ! effective latent heat

      temps = tc
      ps = p
      qs = GEOS_QSAT(temps, ps/100.)
      rvls = qs
 
      if ( qs .ge. qt ) then  ! no condensation
         id_check = 0
         qv = qt
         qc = 0.
         ql = 0.
         qi = 0.
         th = thl !tc/exnerfn(p)
      else                    ! condensation
         do iteration = 1,10
            temps = temps + ( (tc-temps)*cp/leff + qt - rvls )/( cp/leff + ep2*leff*rvls/(r*temps*temps) )
            qs = GEOS_QSAT(temps, ps/100.)
            rvls = qs
         end do
         qc = max( qt-qs, 0._r8 )
         qv = qt - qc
         ql = qc*(1._r8 - nu)
         qi = nu*qc
         th = temps/exnerfn(p)
         if ( abs((temps-(leff/cp)*qc)-tc) .ge. 1._r8 ) then
            id_check = 1
         else
            id_check = 0
         end if
      end if

      return
   end subroutine conden


  subroutine roots(a,b,c,r1,r2,status)
  ! --------------------------------------------------------- !
  ! Subroutine to solve the second order polynomial equation. !
  ! I should check this subroutine later.                     !
  ! --------------------------------------------------------- !
    real, intent(in)  :: a,b,c
    real, intent(out) :: r1,r2
    integer , intent(out) :: status
    real              :: q

    status = 0

    if( a .eq. 0. ) then                            ! Form b*x + c = 0
        if( b .eq. 0. ) then                        ! Failure: c = 0
            status = 1
        else                                           ! b*x + c = 0
            r1 = -c/b
        endif
        r2 = r1
    else
        if( b .eq. 0. ) then                        ! Form a*x**2 + c = 0
            if( a*c .gt. 0. ) then                  ! Failure: x**2 = -c/a < 0
                status = 2  
            else                                       ! x**2 = -c/a 
                r1 = sqrt(-c/a)
            endif
            r2 = -r1
       else                                            ! Form a*x**2 + b*x + c = 0
            if( (b**2 - 4.*a*c) .lt. 0. ) then   ! Failure, no real roots
                 status = 3
            else
                 q  = -0.5*(b + sign(1.0,b)*sqrt(b**2 - 4.*a*c))
                 r1 =  q/a
                 r2 =  c/q
!                 r1 = -0.5*(b + sign(1.0,b)*sqrt(b**2 - 4.*a*c))/a
!                 r2 = -0.5*(b - sign(1.0,b)*sqrt(b**2 - 4.*a*c))/a
            endif
       endif
    endif

    return
  end subroutine roots


   subroutine getbuoy(pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin)
  ! ----------------------------------------------------------- !
  ! Subroutine to calculate integrated CIN [ J/kg = m2/s2 ] and !
  ! 'cinlcl, plfc' if any. Assume 'thv' is linear in each layer !
  ! both for cumulus and environment. Note that this subroutine !
  ! only includes positive CIN in calculation - if there is any !
  ! negative CIN, it is assumed to be zero.    This is slightly !
  ! different from 'single_cin' below, where both positive  and !
  ! negative CIN are included.                                  !
  ! ----------------------------------------------------------- !
    real pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin,frc

    if( thvubot .gt. thv0bot .and. thvutop .gt. thv0top ) then
        plfc = pbot
        return
    elseif( thvubot .le. thv0bot .and. thvutop .le. thv0top ) then 
        cin  = cin - ( (thvubot/thv0bot - 1.) + (thvutop/thv0top - 1.)) * (pbot - ptop) /        &
                     ( pbot/(r*thv0bot*exnerfn(pbot)) + ptop/(r*thv0top*exnerfn(ptop)) )
    elseif( thvubot .gt. thv0bot .and. thvutop .le. thv0top ) then 
        frc  = ( thvutop/thv0top - 1.) / ( (thvutop/thv0top - 1.) - (thvubot/thv0bot - 1.) )
        cin  = cin - ( thvutop/thv0top - 1.) * ( (ptop + frc*(pbot - ptop)) - ptop ) /             &
                     ( pbot/(r*thv0bot*exnerfn(pbot)) + ptop/(r*thv0top*exnerfn(ptop)) )
    else            
        frc  = ( thvubot/thv0bot - 1.) / ( (thvubot/thv0bot - 1.) - (thvutop/thv0top - 1.) )
        plfc = pbot - frc * ( pbot - ptop )
        cin  = cin - ( thvubot/thv0bot - 1.)*(pbot - plfc)/                                         & 
                     ( pbot/(r*thv0bot*exnerfn(pbot)) + ptop/(r*thv0top * exnerfn(ptop)))
    endif

    return
  end subroutine getbuoy

  function single_cin(pbot,thv0bot,ptop,thv0top,thvubot,thvutop)
  ! ------------------------------------------------------- !
  ! Function to calculate a single layer CIN by summing all ! 
  ! positive and negative CIN.                              !
  ! ------------------------------------------------------- ! 
    real :: single_cin
    real :: pbot,thv0bot,ptop,thv0top,thvubot,thvutop 

    single_cin = ( (1.- thvubot/thv0bot) + (1.- thvutop/thv0top)) * ( pbot - ptop ) / &
                 ( pbot/(r*thv0bot*exnerfn(pbot)) + ptop/(r*thv0top*exnerfn(ptop)) )
    return
  end function single_cin


end module uwshcu
