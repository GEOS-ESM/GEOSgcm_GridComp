module uwshcu

!#define UWDIAG 1

  !  use GEOS_Mod, only: write_parallel
  !  use MAPL,     only: MAPL_UNDEF

   use GEOS_UtilsMod, only: GEOS_QSAT, GEOS_DQSAT
   use SHLWPARAMS
   use MAPL_ConstantsMod, only: MAPL_TICE , MAPL_CP   , &
                                MAPL_GRAV , MAPL_ALHS , &
                                MAPL_ALHL , MAPL_ALHF , &
                                MAPL_RGAS , MAPL_H2OMW, &
                                MAPL_AIRMW, MAPL_RVAP , &
                                MAPL_PI, r8 => MAPL_R8

   implicit none

   private
   public compute_uwshcu_inv

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

contains

   real function exnerfn(pressure)
!$acc routine seq
      real, intent(in)              :: pressure   ! in Pa
      exnerfn = (pressure/p00)**rovcp
      return
   end function exnerfn


   subroutine compute_uwshcu_inv(idim, k0, ncnst, dt,pmid0_inv,     & ! INPUT
         zmid0_inv, exnmid0_inv, pifc0_inv, zifc0_inv, exnifc0_inv, &
         dp0_inv, u0_inv, v0_inv, qv0_inv, ql0_inv, qi0_inv,        &
         th0_inv, tke_inv, kpbl_inv, thlsrc_pert,                   & 
         cush, tr0_inv,                                             & ! INOUT
         umf_inv, qvten_inv, qlten_inv, qiten_inv, thten_inv,       & ! OUTPUT
         uten_inv, vten_inv, qrten_inv, qsten_inv, cufrc_inv,       &
         fer_inv, fdr_inv, qldet_inv, qidet_inv, qlsub_inv,         &
         qisub_inv, ndrop_inv, nice_inv,                            & 
#ifdef UWDIAG
         qcu_inv, qlu_inv, qiu_inv, cbmf, qc_inv,                   & ! DIAGNOSTIC ONLY
         cnt_inv, cnb_inv, cin, plcl, plfc, pinv, prel, pbup,       &
         wlcl, qtsrc, thlsrc, thvlsrc, tkeavg, cldtop, wu_inv,      &
         qtu_inv, thlu_inv, thvu_inv, uu_inv, vu_inv, xc_inv,       &
#endif
         dotransport,shlwparams)

      implicit none


      integer, intent(in)   :: idim                     ! Number of columns (IDIM)
      integer, intent(in)   :: k0                       ! Number of levels  (K0) 
      integer, intent(in)   :: ncnst                    ! Number of tracers (ITRCR)
      integer, intent(in)   :: dotransport              ! Transport tracers [1 true] (USE_TRACER_TRANSP_UW)
      real   , intent(in)   :: dt                       ! moist heartbeat [s] (DT_MOIST)

      real,   intent(in)    :: pifc0_inv(idim,k0+1)     !  Environmental pressure at the interfaces [ Pa ] (PLE)
      real,   intent(in)    :: zifc0_inv(idim,k0+1)     !  Environmental height at the interfaces   [ m ]  (ZLE)
      real,   intent(in)    :: exnifc0_inv(idim,k0+1)   !  Exner function at the interfaces (PKE)
      real,   intent(in)    :: pmid0_inv(idim,k0)       !  Environmental pressure at the layer mid-point [ Pa ] (PLO*100.)
      real,   intent(in)    :: zmid0_inv(idim,k0)       !  Environmental height at the layer mid-point [ m ] (ZLO)
      real,   intent(in)    :: exnmid0_inv(idim,k0)     !  Exner function at the layer mid-point (PK)
      real,   intent(in)    :: dp0_inv(idim,k0)         !  Environmental layer pressure thickness [ Pa ] > 0. (DP)
      real,   intent(in)    :: u0_inv(idim,k0)          !  Environmental zonal wind [ m/s ] (U1)
      real,   intent(in)    :: v0_inv(idim,k0)          !  Environmental meridional wind [ m/s ] (V1)
      real,   intent(in)    :: qv0_inv(idim,k0)         !  Environmental water vapor specific humidity [ kg/kg ] (Q1)
      real,   intent(in)    :: ql0_inv(idim,k0)         !  Environmental liquid water specific humidity [ kg/kg ] (QLLS)
      real,   intent(in)    :: qi0_inv(idim,k0)         !  Environmental ice specific humidity [ kg/kg ] (QILS)
      real,   intent(in)    :: th0_inv(idim,k0)         !  Environmental temperature [ K ] (TH1)
      real,   intent(in)    :: tke_inv(idim,k0+1)       !  Turbulent kinetic energy at the interfaces [ m2/s2 ] (TKE)
                                                        !  at the previous time step [ fraction ]
      real, intent(in)    :: kpbl_inv(idim)               !  Height of PBL [ m ] (KPBLSC)
      real, intent(in)    :: thlsrc_pert              !  (THLSRC_PERT)
      real, intent(inout) :: tr0_inv(idim,k0,ncnst)   !  Environmental tracers [ #, kg/kg ] (XHO)
      real, intent(inout) :: cush(idim)               !  Convective scale height [m] (CUSH)

      real, intent(out)   :: umf_inv(idim,k0+1)       !  Updraft mass flux at interfaces [kg/m2/s] (UMF_SC)
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

      type(shlwparam_type), intent(in) :: shlwparams

  !----- Internal variables -----
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
      real              :: tr0(idim,k0,ncnst)       !  Environmental tracers [ #, kg/kg ]
      real              :: umf(idim,0:k0)           !  Updraft mass flux at the interfaces [ kg/m2/s ]
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


!--------- Internal, Diagnostic only ---------
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
      integer :: i



!---------- Indices -----------
      integer           :: k                        !  Vertical index for local fields [ no ] 
      integer           :: k_inv                    !  Vertical index for incoming fields [ no ]
      integer           :: m                        !  Tracer index [ no ]
      integer           :: kpbl(idim)

!!$acc update device(pmid0_inv,zmid0_inv,exnmid0_inv,pifc0_inv,zifc0_inv,exnifc0_inv, &
!!$acc        dp0_inv,u0_inv,v0_inv,qv0_inv,ql0_inv,qi0_inv,th0_inv,tke_inv,kpbl_inv, &
!!$acc        cush,tr0_inv,umf_inv,qvten_inv,qlten_inv,qiten_inv,thten_inv,&
!!$acc        uten_inv,vten_inv,qrten_inv,qsten_inv,cufrc_inv,fer_inv,fdr_inv,qldet_inv,&
!!$acc        qidet_inv,qlsub_inv,qisub_inv,ndrop_inv,nice_inv)

!$acc data create(pmid0,u0,v0,zmid0,exnmid0,dp0,qv0,ql0,qi0,th0,tr0,&
!$acc      tke,pifc0,zifc0,exnifc0,kpbl,umf,qvten,qlten,qiten,sten,&
!$acc      uten,vten,qrten,qsten,cufrc,fer,fdr,qldet,qidet,qlsub,&
!$acc      qisub,ndrop,nice)&
!$acc present(pmid0_inv,zmid0_inv,exnmid0_inv,pifc0_inv,zifc0_inv,exnifc0_inv, &
!$acc dp0_inv,u0_inv,v0_inv,qv0_inv,ql0_inv,qi0_inv,th0_inv,tke_inv,kpbl_inv, &
!$acc cush,tr0_inv,umf_inv,qvten_inv,qlten_inv,qiten_inv,thten_inv,&
!$acc uten_inv,vten_inv,qrten_inv,qsten_inv,cufrc_inv,fer_inv,fdr_inv,qldet_inv,&
!$acc        qidet_inv,qlsub_inv,qisub_inv,ndrop_inv,nice_inv)

      ! flip mid-level variables
!$acc parallel vector_length(128)
!$acc loop gang vector collapse(2) private(k_inv)
      do k = 1, k0
      do i = 1, idim
         k_inv               = k0 + 1 - k
         pmid0(i,k)      = pmid0_inv(i,k_inv)
         u0(i,k)         = u0_inv(i,k_inv)
         v0(i,k)         = v0_inv(i,k_inv)
         zmid0(i,k)      = zmid0_inv(i,k_inv)
         exnmid0(i,k)    = exnmid0_inv(i,k_inv)
         dp0(i,k)        = dp0_inv(i,k_inv)
         qv0(i,k)        = qv0_inv(i,k_inv)
         ql0(i,k)        = ql0_inv(i,k_inv)
         qi0(i,k)        = qi0_inv(i,k_inv)
         th0(i,k)        = th0_inv(i,k_inv)
         do m = 1, ncnst
            tr0(i,k,m)   = tr0_inv(i,k_inv,m)
         enddo
      enddo
      enddo
!$acc end parallel

!$acc parallel vector_length(128)
!$acc loop gang
      ! flip interface variables
      do i = 1, idim
!$acc loop vector
      do k = 0, k0
         tke(i,k) = 0.
         pifc0(i,k) = 0.
         zifc0(i,k) = 0.
         exnifc0(i,k) = 0.
      enddo
!$acc loop vector private(k_inv)
      do k = 0, k0
         k_inv               = k0 - k + 1
         tke(i,k)        = tke_inv(i,k_inv)
         pifc0(i,k)      = pifc0_inv(i,k_inv)
         zifc0(i,k)      = zifc0_inv(i,k_inv)
         exnifc0(i,k)    = exnifc0_inv(i,k_inv)
      end do
      enddo
!$acc end parallel

!$acc parallel vector_length(128)
!$acc loop gang vector
      do i = 1, idim
         kpbl(i) = int(kpbl_inv(i))
      enddo
!$acc end parallel

      call compute_uwshcu( idim,k0, dt, ncnst,pifc0, zifc0, &
           exnifc0, pmid0, zmid0, exnmid0, dp0, u0, v0,     &
           qv0, ql0, qi0, th0, tr0, kpbl, tke, cush, umf,   &
           qvten, qlten, qiten, sten, uten, vten,           &
           qrten, qsten, cufrc, fer, fdr, qldet, qidet,     & 
           qlsub, qisub, ndrop, nice, thlsrc_pert,          &
#ifdef UWDIAG
           qcu, qlu, qiu, cbmf, qc, cnt, cnb,               & ! Diagnostic only
           cin, plcl, plfc, pinv, prel, pbup, wlcl, qtsrc,  &
           thlsrc, thvlsrc, tkeavg, cldtop, wu, qtu,        &
           thlu, thvu, uu, vu, xc, trten,                   & 
#endif
           dotransport,shlwparams )

      ! Reverse again

#ifdef UWDIAG
      cnt_inv(:idim) = k0 + 1 - cnt(:idim)
      cnb_inv(:idim) = k0 + 1 - cnb(:idim)
#endif

!$acc parallel vector_length(128)
!$acc loop gang vector collapse(2) private(k_inv)
      do k = 0, k0
      do i = 1, idim
         k_inv                    = k0 + 1 - k
         umf_inv(i,k_inv)     = umf(i,k)

#ifdef UWDIAG
         wu_inv(i,k_inv)      = wu(i,k)   ! Diagnostic only
         qtu_inv(i,k_inv)     = qtu(i,k)
         thlu_inv(i,k_inv)    = thlu(i,k)
         thvu_inv(i,k_inv)    = thvu(i,k)
         uu_inv(i,k_inv)      = uu(i,k)
         vu_inv(i,k_inv)      = vu(i,k)
#endif
      end do
      enddo
!$acc end parallel

!$acc parallel vector_length(128)
!$acc loop gang vector collapse(2) private(k_inv)
      do k = 1, k0
      do i = 1, idim
         k_inv                    = k0 + 1 - k
         qvten_inv(i,k_inv)   = qvten(i,k)
         qlten_inv(i,k_inv)   = qlten(i,k)
         qiten_inv(i,k_inv)   = qiten(i,k)
         thten_inv(i,k_inv)   = sten(i,k) / (cp*exnmid0(i,k))
         uten_inv(i,k_inv)    = uten(i,k)
         vten_inv(i,k_inv)    = vten(i,k)
         qrten_inv(i,k_inv)   = qrten(i,k)
         qsten_inv(i,k_inv)   = qsten(i,k)
         cufrc_inv(i,k_inv)   = cufrc(i,k)
         fer_inv(i,k_inv)     = fer(i,k)
         fdr_inv(i,k_inv)     = fdr(i,k)
         qldet_inv(i,k_inv)   = qldet(i,k)
         qidet_inv(i,k_inv)   = qidet(i,k)
         qlsub_inv(i,k_inv)   = qlsub(i,k)
         qisub_inv(i,k_inv)   = qisub(i,k)
         ndrop_inv(i,k_inv)   = ndrop(i,k)
         nice_inv(i,k_inv)    = nice(i,k)
#ifdef UWDIAG
         qcu_inv(i,k_inv)     = qcu(i,k)  ! Diagnostic only
         qlu_inv(i,k_inv)     = qlu(i,k)
         qiu_inv(i,k_inv)     = qiu(i,k)
         qc_inv(i,k_inv)      = qc(i,k)
         xc_inv(i,k_inv)      = xc(i,k)
#endif
         if (dotransport.eq.1) then
         do m = 1, ncnst
            tr0_inv(i,k_inv,m)   = tr0(i,k,m)
#ifdef UWDIAG
            trten_inv(i,k_inv,m) = trten(i,k,m)
#endif
         enddo
         endif

      end do
      enddo
!$acc end parallel
!$acc end data
!!!$acc update host(cush,tr0_inv,umf_inv,qvten_inv,qlten_inv,qiten_inv,thten_inv,&
!!!$acc uten_inv,vten_inv,qrten_inv,qsten_inv,cufrc_inv,fer_inv,fdr_inv,qldet_inv,&
!!!$acc qidet_inv,qlsub_inv,qisub_inv,ndrop_inv,nice_inv)

   end subroutine compute_uwshcu_inv


   subroutine compute_uwshcu(idim, k0, dt,ncnst, pifc0_in,zifc0_in,& ! IN
         exnifc0_in, pmid0_in, zmid0_in, exnmid0_in, dp0_in,       &
         u0_in, v0_in, qv0_in, ql0_in, qi0_in, th0_in,             &
         tr0_inout, kpbl_in, tke_in, cush_inout,                   & ! OUT
         umf_out, qvten_out, qlten_out, qiten_out,                 &
         sten_out, uten_out, vten_out, qrten_out,                  &
         qsten_out, cufrc_out, fer_out, fdr_out, qldet_out,        &
         qidet_out, qlsub_out, qisub_out, ndrop_out, nice_out,     &
         thlsrc_pert,                                              &
#ifdef UWDIAG
         qcu_out, qlu_out, qiu_out, cbmf_out, qc_out,              & ! DIAG ONLY
         cnt_out, cnb_out, cinh_out, plcl_out, plfc_out, pinv_out, &
         prel_out, pbup_out, wlcl_out, qtsrc_out, thlsrc_out,      &
         thvlsrc_out, tkeavg_out, cldhgt_out, wu_out, qtu_out,     &
         thlu_out, thvu_out, uu_out, vu_out, xc_out, trten_out,    &
#endif
         dotransport,shlwparams )  

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
      real, intent(in)    :: thlsrc_pert
      integer, intent(in) :: kpbl_in( idim )          ! Boundary layer top layer index

      real, intent(inout) :: cush_inout( idim )       ! Convective scale height [m]
      real, intent(inout) :: tr0_inout(idim,k0,ncnst) !  Environmental tracers [ #, kg/kg ]

      real, intent(out)   :: umf_out(idim,0:k0)       !  Updraft mass flux at the interfaces [ kg/m2/s ]
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

      type(shlwparam_type), intent(in) :: shlwparams

!srf
      !real, intent(in)    :: qtsrc_pert 
      real     :: qtsrc_pert = 0.5e-3
!srf      
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
      real       sigmaw, tkeavg, dpsum, dpi, thvlmin
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

    real, dimension(k0)   :: qv0_s  , ql0_s   , qi0_s   , s0_s    , u0_s    ,          & 
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
    real    :: rpen     ! penetrative entrainment efficiency

    ! ----------------------- !
    ! For lateral entrainment !
    ! ----------------------- !

    real :: rle          !  For critical stopping distance for lateral entrainment [no unit]
    real :: rkm          !  Determine the amount of air that is involved in buoyancy-sorting [no unit]
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
    real :: frc_rasn
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
    rpen             = shlwparams%rpen
    rle       = shlwparams%rle      !  For critical stopping distance for lateral entrainment [no unit]
    rkm       = shlwparams%rkm      !  Determine the amount of air that is involved in buoyancy-sorting [no unit]
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

!!!$acc update device(pifc0_in,zifc0_in,exnifc0_in,pmid0_in,zmid0_in,exnmid0_in,&
!!!$acc        dp0_in,u0_in,v0_in,qv0_in,ql0_in,qi0_in,th0_in,tr0_inout,kpbl_in,&
!!!$acc        tke_in,cush_inout,umf_out,qvten_out,qlten_out,qiten_out,sten_out,&
!!!$acc        uten_out,vten_out,qrten_out,qsten_out,cufrc_out,fer_out,fdr_out,&
!!!$acc        qldet_out,qidet_out,qlsub_out,qisub_out,ndrop_out,nice_out)

!$acc data create(exit_ufrc,exit_wtw,exit_drycore,exit_wu,exit_cufilter,&
!$acc      exit_rei,exit_kinv1,exit_klfck0,exit_klclk0,exit_uwcu,exit_conden,&
!$acc      limit_cinlcl,limit_cin,ind_delcin,limit_rei,limit_shcu,limit_negcon,&
!$acc      limit_ufrc,limit_ppen,limit_emf,limit_cbmf)&
!$acc      present(pifc0_in,zifc0_in,exnifc0_in,pmid0_in,zmid0_in,exnmid0_in,&
!$acc        dp0_in,u0_in,v0_in,qv0_in,ql0_in,qi0_in,th0_in,tr0_inout,kpbl_in,&
!$acc        tke_in,cush_inout,umf_out,qvten_out,qlten_out,qiten_out,sten_out,&
!$acc        uten_out,vten_out,qrten_out,qsten_out,cufrc_out,fer_out,fdr_out,&
!$acc        qldet_out,qidet_out,qlsub_out,qisub_out,ndrop_out,nice_out)
    ! ------------------------------------------------------- !
    ! Initialize output variables defined for all grid points !
    ! ------------------------------------------------------- !

!$acc parallel vector_length(128)
!$acc loop gang 
    do i = 1, idim
!$acc loop vector
    do k = 1, k0
    umf_out(i,k)          = 0.0
    cufrc_out(i,k)         = 0.0
    fer_out(i,k)           = 0.0
    fdr_out(i,k)           = 0.0
    qldet_out(i,k)         = 0.0
    qidet_out(i,k)         = 0.0
    qlsub_out(i,k)         = 0.0
    qisub_out(i,k)         = 0.0
    ndrop_out(i,k)         = 0.0
    nice_out(i,k)         = 0.0
    enddo
    umf_out(i,0)          = 0.0

#ifdef UWDIAG
    cbmf_out(:idim)              = 0.0
    cinh_out(:idim)              = -987
    cinlclh_out(:idim)           = -987
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

    wu_out(:idim,0:k0)          = -987
    qtu_out(:idim,0:k0)         = -987
    thlu_out(:idim,0:k0)        = -987
    thvu_out(:idim,0:k0)        = -987
    uu_out(:idim,0:k0)          = -987
    vu_out(:idim,0:k0)          = -987
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

    exit_UWCu(i)             = 0.0
    exit_conden(i)           = 0.0
    exit_klclk0(i)           = 0.0
    exit_klfck0(i)           = 0.0
    exit_ufrc(i)             = 0.0
    exit_wtw(i)              = 0.0
    exit_drycore(i)          = 0.0
    exit_wu(i)               = 0.0
    exit_cufilter(i)         = 0.0
    exit_kinv1(i)            = 0.0
    exit_rei(i)              = 0.0

    limit_shcu(i)            = 0.0
    limit_negcon(i)          = 0.0
    limit_ufrc(i)            = 0.0
    limit_ppen(i)            = 0.0
    limit_emf(i)             = 0.0
    limit_cinlcl(i)          = 0.0
    limit_cin(i)             = 0.0
    limit_cbmf(i)            = 0.0
    limit_rei(i)             = 0.0

    ind_delcin(i)            = 0.0
    enddo
!$acc end parallel


      !========================
      !   Start column loop
      !========================
!!!$omp private(id_exit,pifc0,zifc0,pmid0,zmid0,dp0,u0,v0,qv0,ql0,qi0,tke,cush, &
!!!$omp tr0,exnmid0,exnifc0,t0,s0,qt0,thl0,thvl0,ssthl0,ssqt0,ssu0,ssv0,sstr0,  &
!!!$omp thl0top,thl0bot,qt0bot,qt0top,thvubot,thvutop,thj,qvj,qlj,qij,thvj,tj,  &
!!!$omp thv0j,rhomid0j,rhoifc0j,qse,id_check,thv0bot,thvl0bot,thv0top, &
!!!$omp thvl0top,qv0_o,ql0_o,qi0_o,t0_o,s0_o,u0_o,v0_o,qt0_o,thl0_o,thvl0_o, &
!!!$omp ssthl0_o,ssqt0_o,ssu0_o,ssv0_o,thv0bot_o,thv0top_o,thvl0bot_o,thvl0top_o, &
!!!$omp tr0_o,trten_o,sstr0_o,umf,emf,slflx,qtflx,uflx,vflx,qvten,qlten,sten,uten, &
!!!$omp vten,qrten,qsten,dwten,diten,cufrc,qcu,qlu,qiu,fer,xco,cin,cinlcl,cbmf, &
!!!$omp wcrit,winv,wlcl,ufrcinv,ufrclcl,exql,exqi,ppen, &
!!!$omp thl0lcl,qt0lcl,thv0lcl,thv0rel,rho0inv,autodet,thlu_top,qtu_top,qlu_top, &
!!!$omp qiu_top,qlu_mid,qiu_mid,exntop,aquad,bquad,cquad,xc1,xc2,excessu,excess0, &
!!!$omp xsat,xs1,xs2,bogbot,bogtop,delbog,drage,expfac,rcwp,rlwp,riwp,qcubelow, &
!!!$omp qlubelow,qiubelow,qs,qsat_arg,qsat_pe,pe,xsrc,xmean,xtop,xbot,xflx,qiten,fdr, &
!!!$omp qc,qc_l,qc_i,cnt,cnb,qtten,slten,ufrc,thlu,qtu,uu,vu,wu,thvu,thlu_emf, &
!!!$omp qtu_emf,uu_emf,ufrcinvbase,winvbase,emfkbup,cbmflimit,uemf,comsub,qlten_sink, &
!!!$omp qiten_sink,qlten_det,qiten_det,vu_emf,trflx,trten,tru,tru_emf,  &
!!!$omp cldhgt,scaleh,tscaleh,cridis,sigmaw,tkeavg,dpsum,dpi,thvlmin,kinv,thlsrc, &
!!!$omp qtsrc,usrc,vsrc,thvlsrc,plcl,plfc,prel,wrel,klfc,cin_i,cin_f,del_CIN, &
!!!$omp ke,alpha,kinv_o,klcl_o,klfc_o,klcl,usrc_o, &
!!!$omp vsrc_o,plcl_o,plfc_o,thv0lcl_o,trsrc_o,cinlcl_i,cinlcl_f,del_cinlcl, &
!!!$omp qlten_sub,qiten_sub,umf_s ,slflx_s,qtflx_s,ufrc_s ,uflx_s,vflx_s,qvten_s, &
!!!$omp qlten_s,qiten_s,qrten_s,qsten_s,sten_s,uten_s,vten_s,cufrc_s,qcu_s,qlu_s, &
!!!$omp qiu_s,fer_s,fdr_s,xc_s,qc_s,qtten_s,slten_s,qldet_s,qidet_s,qlsub_s,qisub_s, &
!!!$omp tkeavg_o,thvlmin_o,qtsrc_o,thvlsrc_o,thlsrc_o,cush_s,cin_s,cbmf_s,cnt_s, &
!!!$omp cnb_s,krel,mu,mumin0,mumin2,mulcl,mulclstar,ee2,ud2,wtw,wtwb,uplus,vplus, &
!!!$omp dpe,exne,thvebot,thle,qte,ue,ve,thlue,qtue,wue,trsrc,tre,iter_scaleh, &
!!!$omp iter_xc,kbup,kpen,kk,k,i,kp1,km1,mm,m,xc,thlxsat,qtxsat,thvxsat,x_cu, &
!!!$omp x_en,thv_x0,thv_x1,rei,forcedCu,thlten_sub,qtten_sub,nlten_sub,niten_sub, &
!!!$omp thl_prog,qt_prog,uf,vf,qc_lm,nc_lm,nc_im,ql_emf_kbup,status,qc_im, &
!!!$omp qi_emf_kbup,nl_emf_kbup,ni_emf_kbup,totsink,qv0_star,ql0_star,qi0_star, &
!!!$omp s0_star,trmin,trflx_d,trflx_u,pdelx,dum,qv0_s,ql0_s,qi0_s,s0_s,u0_s, &
!!!$omp v0_s,t0_s,tr0_s,cinlcl_s) &
!#ifdef UWDIAG
!#endif

!$omp parallel do default(private)   &
!$omp shared(iter_cin,niter_xc,use_CINcin,use_self_detrain,use_momenflx, &
!$omp use_cumpenent,scverbose,rpen,rle,rkm,rkfre,rmaxfrac,mumin1,rbuoy,rdrag,  &
!$omp epsvarw,PGFc,criqc,frc_rasn,rdrop,ixcldice,ixcldliq,ixnumice,ixnumliq, &
!$omp pifc0_in,zifc0_in,exnifc0_in,pmid0_in,zmid0_in,exnmid0_in,       &
!$omp dp0_in,u0_in,v0_in,qv0_in,ql0_in,qi0_in,th0_in,tr0_inout,tke_in,        &
!$omp kpbl_in,cush_inout,umf_out,qvten_out,qlten_out,qiten_out,   &
!$omp sten_out,uten_out,vten_out,qrten_out,qsten_out, &
!$omp cufrc_out,qldet_out,qidet_out,qlsub_out,qisub_out,ndrop_out,      &
!$omp nice_out,fer_out,fdr_out,    &
!$omp exit_ufrc,exit_wtw,exit_wu,exit_cufilter,limit_rei,  &
!$omp exit_kinv1,exit_klfck0,exit_klclk0,exit_uwcu,exit_conden,exit_rei,      &
!$omp limit_cinlcl,limit_cin,ind_delcin,limit_shcu,dt,thlsrc_pert,limit_negcon,&
#ifdef UWDIAG
!$omp cbmf_out,cinh_out,cinlclh_out,cldhgt_out,qcu_out,qlu_out,qiu_out,      &
!$omp qc_out,cnt_out,cnb_out,xc_out,ufrcinvbase_out,ufrclcl_out,           &
!$omp winvbase_out,wlcl_out,plcl_out,pinv_out,plfc_out,prel_out,pbup_out,    &
!$omp ppen_out,qtsrc_out,thlsrc_out,thvlsrc_out,emfkbup_out,      &
!$omp cbmflimit_out,zinv_out,rcwp_out,rlwp_out,emfkbup_out,cbmflimit_out, &
!$omp zinv_out,rcwp_out,rlwp_out,riwp_out,tkeavg_out,wu_out,qtu_out        &
!$omp thlu_out,thvu_out,uu_out,vu_out,qtu_emf_out,thlu_emf_out,uu_emf_out, &
!$omp vu_emf_out,uemf_out,dwten_out,diten_out,trten_out,trflx_out,tru_out,    &
!$omp tru_emf_out,excessu_arr_out,excess0_arr_out,xc_arr_out,aquad_arr_out, &
!$omp bquad_arr_out,cquad_arr_out,bogbot_arr_out,bogtop_arr_out
#endif
!$omp  limit_ufrc,limit_ppen,limit_emf,limit_cbmf,exit_drycore,idim,k0,ncnst, &
!$omp dotransport)

!$acc parallel vector_length(32)
!$acc loop gang
      do i = 1, idim

         id_exit = .false.
!$acc loop vector
         do k = 1,k0
         pifc0(k)     = pifc0_in(i,k)
         zifc0(k)     = zifc0_in(i,k)
         pmid0(k)      = pmid0_in(i,k)
         zmid0(k)      = zmid0_in(i,k)
         dp0(k)        = dp0_in(i,k)
         u0(k)         = u0_in(i,k)
         v0(k)         = v0_in(i,k)
         qv0(k)        = qv0_in(i,k)
         ql0(k)        = ql0_in(i,k)
         qi0(k)        = qi0_in(i,k)
         tke(k)       = tke_in(i,k)
         enddo
         pifc0(0)     = pifc0_in(i,0)
         zifc0(0)     = zifc0_in(i,0)
!         pblh            = pblh_in(i)
         cush            = cush_inout(i)

         if (dotransport.eq.1) then
!$acc loop vector collapse(2)
         do m = 1,ncnst   ! loop over tracers
         do k = 1,k0
            tr0(k,m) = tr0_inout(i,k,m)
         end do
         enddo
         endif

         !------------------------------------------------------!
         ! Compute basic thermodynamic variables directly from  !
         ! input variables for each column                      !
         !------------------------------------------------------!

         ! Compute internal environmental variables
!$acc loop vector
         do k = 1,k0
         exnmid0(k) = exnmid0_in(i,k)
         exnifc0(k) = exnifc0_in(i,k)
         t0(k)      = th0_in(i,k) * exnmid0(k)
         s0(k)      = g*zmid0(k) + cp*t0(k)
         qt0(k)     = qv0(k) + ql0(k) + qi0(k)
         thl0(k)    = ( t0(k) - xlv*ql0(k)/cp - xls*qi0(k)/cp ) / exnmid0(k)
         thvl0(k)   = ( 1. + zvir*qt0(k) )*thl0(k)
         enddo
         exnifc0(0) = exnifc0_in(i,0)

         ! Compute slopes of environmental variables in each layer

         ssthl0    = slope( k0, thl0, pmid0 )
         ssqt0     = slope( k0, qt0 , pmid0 )
         ssu0      = slope( k0, u0  , pmid0 )
         ssv0      = slope( k0, v0  , pmid0 )
         if (dotransport.eq.1) then
!!!         do m = 1,ncnst    ! loop over tracers
!!!            sstr0(:k0,m) = slope1( k0, tr0(:k0,m), pmid0 )
            call slope1(k0, ncnst, tr0, pmid0, sstr0)
!!!         end do
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
                  !call write_parallel('------- UW ShCu: Exit, conden')
               end if
!!!               write(*,*) "Print1"
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
                   !call write_parallel('------- UW ShCu: Exit, conden')
                 end if
!!!                 write(*,*) "Print2"
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

!!!         write(*,*) "Print3" 
         ! ---------------------------------------------- !
         ! Initialize output variables at each grid point !
         ! ---------------------------------------------- !

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
!      qldet(:k0)         = 0.0
!      qidet(:k0)         = 0.0
      qc_l(:k0)          = 0.0
      qc_i(:k0)          = 0.0
      qtten(:k0)         = 0.0
      slten(:k0)         = 0.0   
      comsub(:k0)        = 0.0
      qlten_sink(:k0)    = 0.0
      qiten_sink(:k0)    = 0.0
      qlten_det(:k0)     = 0.0
      qiten_det(:k0)     = 0.0 

      ufrc(0:k0)         = 0.0  
      umf(0:k0)          = 0.0
      emf(0:k0)          = 0.0
      slflx(0:k0)        = 0.0
      qtflx(0:k0)        = 0.0
      uflx(0:k0)         = 0.0
      vflx(0:k0)         = 0.0

      thlu(0:k0)         = -987
      qtu(0:k0)          = -987
      uu(0:k0)           = -987
      vu(0:k0)           = -987
      wu(0:k0)           = -987
      thvu(0:k0)         = -987
      thlu_emf(0:k0)     = -987
      qtu_emf(0:k0)      = -987
      uu_emf(0:k0)       = -987
      vu_emf(0:k0)       = -987
      uemf(0:k0)         = 0.0
 
      cnt                = real(k0)
      cnb                = 0.0
      cin                = 0.0
      cinlcl             = 0.0
      cbmf               = 0.0     
      ufrcinvbase         = 0.0
      ufrclcl             = 0.0
      winvbase            = 0.0
      wlcl                = 0.0
      emfkbup             = 0.0 
      cbmflimit           = 0.0

!      nlten_sink(:k0)    = 0.0
!      niten_sink(:k0)    = 0.0 

      if (dotransport.eq.1) then
!$acc loop vector collapse(2)
      do m = 1, ncnst
      do k = 0,k0
         trflx(k,m)   = 0.0
         tru(k,m)     = 0.0
         tru_emf(k,m) = 0.0
      enddo
      enddo
!$acc loop vector collapse(2)
      do m = 1, ncnst
      do k = 1,k0
         trten(k,m)    = 0.0
      enddo
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

!           do k = k0-1, 1, -1 
!              if( (pblh + 5. - zifc0(k))*(pblh + 5. - zifc0(k+1)) .lt. 0. ) then
!                  kinv = k + 1   ! should this be k?
!                  go to 15
!              endif 
!           end do
!           kinv = 1

           ! invert kpbl index
           if (kpbl_in(i).gt.k0/2) then
             kinv = k0 - kpbl_in(i) + 1
           else
             kinv = 5
           end if
           !print *, kpbl_in(i), kinv, k0, k0/4, k0/2
!           print *,'kinv=',kinv

15         continue    

           if( kinv .le. 1 ) then        
              exit_kinv1(i) = 1.
              id_exit = .true.
              if (scverbose) then
                !call write_parallel('------- UW ShCu: Exit, kinv<=1')
              end if
!!!              write(*,*) "Print4"
              go to 333
           endif

       end do  ! iter loop

!!!       write(*,*) "Print7"
333    if (id_exit) then

!!!       write(*,*) "Print8"
         exit_uwcu(i) = 1.
         if (scverbose) then
           !call write_parallel('------- UW ShCu: Exited!')
         end if

     ! --------------------------------------------------------------------- !
     ! Initialize output variables when cumulus convection was not performed.!
     ! --------------------------------------------------------------------- !
!$acc loop vector
     do k = 1,k0     
     umf_out(i,k)             = 0.   
     qvten_out(i,k)            = 0.
     qlten_out(i,k)            = 0.
     qiten_out(i,k)            = 0.
     sten_out(i,k)             = 0.
     uten_out(i,k)             = 0.
     vten_out(i,k)             = 0.
     qrten_out(i,k)            = 0.
     qsten_out(i,k)            = 0.
     cufrc_out(i,k)            = 0.
     qldet_out(i,k)            = 0.
     qidet_out(i,k)            = 0.

     fer_out(i,k)             = 0.
     fdr_out(i,k)             = 0.
     enddo
     umf_out(i,0)             = 0.
     cush_inout(i)               = -1.

       end if

     end do ! column i loop
!$acc end parallel
!$acc end data
!!!$acc update host(tr0_inout,cush_inout,umf_out,qvten_out,qlten_out,qiten_out,sten_out,&
!!!$acc        uten_out,vten_out,qrten_out,qsten_out,cufrc_out,fer_out,fdr_out,&
!!!$acc        qldet_out,qidet_out,qlsub_out,qisub_out,ndrop_out,nice_out)

     return

   end subroutine compute_uwshcu


   function slope(k0,field,p0)
!$acc routine seq
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

   subroutine slope1(k0,ncnst,field,p0,sstr0)
!$acc routine vector
      integer, intent(in):: k0,ncnst
      real, intent(out):: sstr0(k0,ncnst)
      real, intent(in) :: field(k0,ncnst)
      real, intent(in) :: p0(k0)
      real             :: below
      real             :: above
      integer            :: k,m

         below = ( field(2,m) - field(1,m) ) / ( p0(2) - p0(1) )
!$acc loop vector private(below,above)
      do m = 1, ncnst
!$acc loop seq
      do k = 2, k0
         above = ( field(k,m) - field(k-1,m) ) / ( p0(k) - p0(k-1) )  
         if ( above .gt. 0. ) then
            sstr0(k-1,m) = max(0.,min(above,below))
         else
            sstr0(k-1,m) = min(0.,max(above,below))
         end if
         below = above
      end do
      sstr0(k0,m) = sstr0(k0-1,m)
      enddo

      return
   end subroutine slope1



  function qsinvert(qt,thl,ps_in)
!$acc routine seq
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
          !call write_parallel('Source air is too dry and pLCL is set to psmin in uwshcu.F90')
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
       nu       =  max(min((268._r8 - Ts)/20._r8,1.0_r8),0.0_r8)        
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
!$acc routine seq
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
!$acc routine seq
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
!$acc routine seq
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
!$acc routine seq
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
!$acc routine seq
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
!$acc routine seq
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
      integer            :: iteration
      real*8               :: qs    ! Saturation specific humidity

      tc = thl*exnerfn(p)

      nu = max(min((268._r8-tc)/20._r8,1.0_r8),0.0_r8)    ! ice fraction of condensate
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
!$acc routine seq
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
!$acc routine seq
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
!$acc routine seq
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
