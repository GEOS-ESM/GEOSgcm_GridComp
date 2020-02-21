module edmf_mod

use MAPL_ConstantsMod, only: mapl_grav, mapl_cp, mapl_alhl, mapl_p00, mapl_epsilon, mapl_alhs, mapl_kappa
use MAPL_Mod, only:          mapl_undef

implicit none

real, parameter ::     &
     Wa       = 1.,    &
     Wb       = 1.5,   &
     WSTARmin = 1.e-3, &
     zpblmin  = 100.,  &
     onethird = 1./3.

contains

!
! A_star_closure
!
subroutine A_star_closure(IM, JM, LM, zle, zlo, thv, & ! in
                          izsl, A_star)                ! out

  integer, intent(in)                     :: IM, JM, LM
  real, dimension(IM,JM,LM), intent(in)   :: zlo, thv
  real, dimension(IM,JM,0:LM), intent(in) :: zle
  integer, dimension(IM,JM), intent(out)  :: izsl
  real, dimension(IM,JM,LM), intent(out)  :: A_star

  integer :: i, j, k, km1, kp1
  real :: dthvdz, foo
  real, dimension(IM,JM) :: A_star_sum
  logical, dimension(IM,JM) :: conv_flag, sl_flag

  A_star(:,:,:)   = 0.
  A_star_sum(:,:) = 0.

  ! Determine if column is convective (surface-driven) and compute A_star
  ! (before normalization) at first model level
  do j = 1,JM
  do i = 1,IM
     conv_flag(i,j) = thv(i,j,LM-1) - thv(i,j,LM) < 0.
     sl_flag(i,j)   = conv_flag(i,j)
  end do
  end do

  ! Compute A_star (before normalization) and find top of surface layer
  do k = LM,1,-1
     km1 = k - 1
     kp1 = k + 1

     do j = 1,JM
     do i = 1,IM
        if ( sl_flag(i,j) ) then
           dthvdz = ( thv(i,j,km1) - thv(i,j,k) )/( zlo(i,j,km1) - zlo(i,j,k) )

           if ( dthvdz < 0. ) then
              A_star(i,j,k)   = -sqrt(zle(i,j,km1))*dthvdz
              A_star_sum(i,j) = A_star_sum(i,j) + A_star(i,j,k)*( zle(i,j,km1) - zle(i,j,k) )
           else
              izsl(i,j) = k

              sl_flag(i,j) = .false.
           end if
        end if
     end do
     end do
  end do

  ! Normalize A_star
  ! Note: this loop needs to be refactored for column-major order
  do j = 1,JM
  do i = 1,IM
     if ( conv_flag(i,j) ) then
        do k = LM,izsl(i,j)+1,-1
           A_star(i,j,k) = A_star(i,j,k)/A_star_sum(i,j)
        end do
     end if
  end do
  end do

  ! Test
  do j = 1,JM
  do i = 1,IM
     foo = 0.
     if ( conv_flag(i,j) ) then
        do k = LM,izsl(i,j),-1
           foo = foo + A_star(i,j,k)*( zle(i,j,km1) - zle(i,j,k) ) 
           write(*,*) k, foo, A_star(i,j,k), zle(i,j,km1) - zle(i,j,k)
        end do
     end if
  end do
  end do

  

end subroutine A_star_closure

!
! edmf
!
subroutine run_edmf(IM, JM, LM, numup, iras, jras, &                                ! in
                    edmf_discrete_type, edmf_implicit, thermal_plume_flag, &        ! in
                    th00, dt, z, zle, ple, rhoe, exf, &                             ! in
                    u, v, thl, thv, qt, qv, ql, qi, &                               ! in         
                    ustar, sh, evap, ice_ramp, &                                    ! in
                    pwmin, pwmax, AlphaW, AlphaQT, AlphaTH, &                       ! in
                    ET, L0, ENT0, EDfac, EntWFac, &                                 ! in
                    zpbl, &                                                         ! inout
                    edmfdrya, edmfmoista, &                                         ! out
                    edmfdryw, edmfmoistw, &                                         ! out
                    edmfdryqt, edmfmoistqt, &                                       ! out
                    edmfdrythl, edmfmoistthl, &                                     ! out
                    edmfdryu, edmfmoistu,  &                                        ! out
                    edmfdryv, edmfmoistv,  &                                        ! out
                    edmfmoistqc, &                                                  ! out
                    ae, awu, awv, aw, aws, awqv, awql, awqi, &                      ! out (for solver)
                    whl_mf, wqt_mf, wthv_mf, &                                      ! out (for MYNN-EDMF)
                    buoyf, mfw2, mfw3, mfqt3, mfwqt, mfqt2, mfhl2, mfqthl, mfwhl, & ! out (for SHOC)
                    au, wu, Mu, E, D, hle, qte)                                     ! out
  
  ! Inputs
  integer, intent(in)                     :: IM, JM, LM, numup, edmf_discrete_type, edmf_implicit, thermal_plume_flag, ET
  integer, dimension(IM,JM), intent(in)   :: iras, jras
  real, dimension(IM,JM,LM), intent(in)   :: u, v, thl, qt, thv, qv, ql, qi, z, exf
  real, dimension(IM,JM,0:LM), intent(in) :: zle, ple, rhoe
  real, dimension(IM,JM), intent(in)      :: ustar, sh, evap, L0
  real, dimension(IM,JM), intent(inout)   :: zpbl
  real, intent(in)                        :: th00, ice_ramp, dt, EntWFac, ENT0, pwmin, pwmax, AlphaW, AlphaQT, AlphaTH, EDfac
 
  ! Outputs
  real, dimension(IM,JM,0:LM), intent(out) :: edmfdrya, edmfmoista, edmfdryw, edmfmoistw, &
                                              edmfdryqt, edmfmoistqt, edmfdrythl, edmfmoistthl, &
                                              edmfdryu, edmfmoistu, edmfdryv, edmfmoistv, &
                                              edmfmoistqc, &
                                              ae, aw, aws, awqv, awql, awqi, awu, awv, &
                                              whl_mf, wqt_mf, wthv_mf, au, Mu, wu
  real, dimension(IM,JM,LM), intent(out)   :: buoyf, mfw2, mfw3, mfqt3, mfqt2, mfwqt, &
                                              mfhl2, mfqthl, mfwhl, E, D, hle, qte

  real, dimension(numup,IM,JM) :: upw, upthl, upqt, upql, upqi, upa, upu, upv, upthv

  integer :: i, j, k, km1, kp1, iup

  real :: wthv, wstar, qstar, thstar, sigmaw, sigmaqt, sigmath, z0, wmin, wmax, wlv, wtv, wp, &
          B, QTn, THLn, THVn, QCn, Un, Vn, Wn2, EntEXP, EntEXPU, EntW, wf, &
          stmp, ltm, QTsrfF, THVsrfF, mft, mfthvt, dzle, ifac, test, mft_work, mfthvt_work, &
          au_full, thlu_full, qtu_full

  ! Temporary (too slow; need to figure out how random number generator works)
  integer, dimension(numup,IM,JM,LM)  :: enti
  real, dimension(numup,IM,JM,LM)     :: entf, ent

  real, dimension(IM,JM,LM) :: thlu, qtu

  real, dimension(IM,JM,0:LM) :: aw2, ahl2, aqt2, aw3, aqt3, aqthl, &
                                 ui, vi, thvi, qvi, qli, qii, exfh, thli, qti 

  integer, dimension(2)  :: seedmf, the_seed

  ! New stuff
  logical, dimension(numup,IM,JM) :: active_updraft

  !
  ! Initialize arrays
  !

  ! Updraft statistics
  edmfdrya     = 0.
  edmfmoista   = 0.
  edmfdryw     = 0.
  edmfmoistw   = 0.
  edmfdryqt    = 0.
  edmfmoistqt  = 0.
  edmfdrythl   = 0.
  edmfmoistthl = 0.
  edmfdryu     = 0.
  edmfmoistu   = 0.
  edmfdryv     = 0.
  edmfmoistv   = 0.
  edmfmoistqc  = 0.

  ! Outputs needed for solver 
  aw    = 0.
  awu   = 0.
  awv   = 0.
  aws   = 0.
  awqv  = 0.
  awql  = 0.
  awqi  = 0.

  ! Outputs needed for MYNN-EDMF
  whl_mf  = 0.
  wqt_mf  = 0.
  wthv_mf = 0.

  ! 
  aw2    = 0.
  ahl2   = 0.
  aqt2   = 0.
  aqthl = 0.
  aw3    = 0.
  aqt3   = 0.

  ! Outputs needed for Moist
  buoyf  = 0.
  mfw2   = 0.
  mfw3   = 0.
  mfqt3  = 0.
  mfqt2  = 0.
  mfwqt  = 0.
  mfhl2  = 0.
  mfqthl = 0.
  mfwhl  = 0.

  ae = EDfac 

  au = 0.
  wu = 0.
  Mu = 0.
  E  = 0.
  D  = 0.

  ! Loop across surface tiles
  do j = 1,JM
  do i = 1,IM
     seedmf(1) = 1000000*( 100*thl(i,j,1) - INT(100*thl(i,j,1)) )
     seedmf(2) = 1000000*( 100*thl(i,j,2) - INT(100*thl(i,j,2)) )

     the_seed(1) = seedmf(1)*iras(i,j) + seedmf(2)*jras(i,j)
     the_seed(2) = seedmf(1)*jras(i,j) + seedmf(2)*iras(i,j)
     the_seed(1) = the_seed(1)*seedmf(1)/( seedmf(2) + 10 )
     the_seed(2) = the_seed(2)*seedmf(1)/( seedmf(2) + 10 )
     if(the_seed(1) == 0) the_seed(1) =  5
     if(the_seed(2) == 0) the_seed(2) = -5

     do k = 1,LM
        km1 = k - 1
        do iup = 1,numup
           entf(iup,i,j,k) = ( zle(i,j,km1) - zle(i,j,k) )/L0(i,j)
        end do
     end do

     ! Treat TOA and surface properties the same as the nearest grid cell
     ui(i,j,0)   = u(i,j,1)
     vi(i,j,0)   = v(i,j,1)
     thli(i,j,0) = thl(i,j,1)
     qti(i,j,0)  = qt(i,j,1)
     qvi(i,j,0)  = qv(i,j,1)
     qli(i,j,0)  = ql(i,j,1)
     qii(i,j,0)  = qi(i,j,1)
     thvi(i,j,0) = thv(i,j,1)

     ui(i,j,LM)   = u(i,j,LM)
     vi(i,j,LM)   = v(i,j,LM)
     thli(i,j,LM) = thl(i,j,LM)
     qti(i,j,LM)  = qt(i,j,LM)
     qvi(i,j,LM)  = qv(i,j,LM)
     qli(i,j,LM)  = ql(i,j,LM)
     qii(i,j,LM)  = qi(i,j,LM)
     thvi(i,j,LM) = thv(i,j,LM)

     !
     exfh(i,j,LM) = (ple(i,j,LM)/mapl_p00)**mapl_kappa
  end do
  end do

  ! Interpolate to half levels
  do k = 1,LM-1
     kp1 = k + 1

     do j = 1,JM
     do i = 1,IM
        exfh(i,j,k) = (ple(i,j,k)/mapl_p00)**mapl_kappa

        if (edmf_discrete_type == 0) then
           if (edmf_implicit == 0) then
              ifac = ( zle(i,j,k) - z(i,j,kp1) )/( z(i,j,k) - z(i,j,kp1) )

              ui(i,j,k)   = u(i,j,kp1)   + ifac*( u(i,j,k)   - u(i,j,kp1) )
              vi(i,j,k)   = v(i,j,kp1)   + ifac*( v(i,j,k)   - v(i,j,kp1) )
              thli(i,j,k) = thl(i,j,kp1) + ifac*( thl(i,j,k) - thl(i,j,kp1) )
              qti(i,j,k)  = qt(i,j,kp1)  + ifac*( qt(i,j,k)  - qt(i,j,kp1) )
              qvi(i,j,k)  = qv(i,j,kp1)  + ifac*( qv(i,j,k)  - qv(i,j,kp1) )
              qli(i,j,k)  = ql(i,j,kp1)  + ifac*( ql(i,j,k)  - ql(i,j,kp1) )
              qii(i,j,k)  = qi(i,j,kp1)  + ifac*( qi(i,j,k)  - qi(i,j,kp1) )
              thvi(i,j,k) = thv(i,j,kp1) + ifac*( thv(i,j,k) - thv(i,j,kp1) )
           else ! This is temporary until I fix the interplation for implicit mass flux discretization
              ui(i,j,k)   = 0.5*( u(i,j,kp1)   + u(i,j,k) )
              vi(i,j,k)   = 0.5*( v(i,j,kp1)   + v(i,j,k) )
              thli(i,j,k) = 0.5*( thl(i,j,kp1) + thl(i,j,k) )
              qti(i,j,k)  = 0.5*( qt(i,j,kp1)  + qt(i,j,k) )
              qvi(i,j,k)  = 0.5*( qv(i,j,kp1)  + qv(i,j,k) )
              qli(i,j,k)  = 0.5*( ql(i,j,kp1)  + ql(i,j,k) )
              qii(i,j,k)  = 0.5*( qi(i,j,kp1)  + qi(i,j,k) )
              thvi(i,j,k) = 0.5*( thv(i,j,kp1) + thv(i,j,k) )
           end if
        else
           ui(i,j,k)   = u(i,j,k)
           vi(i,j,k)   = v(i,j,k)
           thli(i,j,k) = thl(i,j,k)
           qti(i,j,k)  = qt(i,j,k)
           qvi(i,j,k)  = qv(i,j,k)
           qli(i,j,k)  = ql(i,j,k)
           qii(i,j,k)  = qi(i,j,k)
           thvi(i,j,k) = thv(i,j,k)
        end if
     end do
     end do
  end do

  ! Initialize updrafts
  do j = 1,JM
  do i = 1,IM
     ent = 0.

     zpbl(i,j) = max( zpbl(i,j), zpblmin )
     wthv      = sh(i,j)/mapl_cp + mapl_epsilon*thv(i,j,LM)*evap(i,j)

     if (wthv > 0.) then
        wstar  = max( wstarmin, (mapl_grav/th00*wthv*zpbl(i,j))**onethird )
        qstar  = evap(i,j)/wstar
        thstar = wthv/wstar

        sigmaW  = AlphaW*wstar
        sigmaQT = AlphaQT*max( qstar, 0. )
        sigmaTH = AlphaTH*max( thstar, 0. )

        wmin = sigmaW*pwmin
        wmax = sigmaW*pwmax

        do iup = 1,numup
           active_updraft(iup,i,j) = .true.

           upw(iup,i,j)   = 0.
           upthl(iup,i,j) = 0.
           upthv(iup,i,j) = 0.
           upqt(iup,i,j)  = 0.
           upa(iup,i,j)   = 0.
           upu(iup,i,j)   = 0.
           upv(iup,i,j)   = 0.
           upqi(iup,i,j)  = 0.
           upql(iup,i,j)  = 0.

           wlv = wmin + ( wmax - wmin )/(real(numup))*( real(iup) - 1. )
           wtv = wmin + ( wmax - wmin )/(real(numup))*real(iup)
       
!           upw(iup,i,j) = min( 0.5*( wlv + wtv ), 5. )  ! npa
           upw(iup,i,j) = 0.5*( wlv + wtv ) 
           upa(iup,i,j) = 0.5*erf(wtv/(sqrt(2.)*sigmaW)) - 0.5*erf(wlv/(sqrt(2.)*sigmaW))
       
           upu(iup,i,j) = u(i,j,LM)
           upv(iup,i,j) = v(i,j,LM)
       
           upqt(iup,i,j)  = qt(i,j,LM)  + 0.32*upw(iup,i,j)*sigmaQT/sigmaW
           upthv(iup,i,j) = thv(i,j,LM) + 0.58*upw(iup,i,j)*sigmaTH/sigmaW

           ! For stability make sure that the surface mass-fluxes are not more than their 
           ! values computed from the surface scheme 
           QTsrfF  = QTsrfF  + upa(iup,i,j)*upw(iup,i,j)*( upqt(iup,i,j)  - qt(i,j,LM) )
           THVsrfF = THVsrfF + upa(iup,i,j)*upw(iup,i,j)*( upthv(iup,i,j) - thv(i,j,LM) )

           ! Change surface THV so that the fluxes from the mass flux equal prescribed values
           if ( THVsrfF > wthv ) then
              upthv(iup,i,j) = ( upthv(iup,i,j) - thv(i,j,LM) )*wthv/THVsrfF + thv(i,j,LM)
           end if

           ! Change surface QT so that the fluxes from the mass flux equal prescribed values 
           ! We do not need to worry about the negative values as they should not exist
           if ( ( QTsrfF > evap(i,j) ) .and. ( evap(i,j) > 0. ) )  then
              upqt(iup,i,j) = ( upqt(iup,i,j) - qt(i,j,LM) )*wthv/QTsrfF + qt(i,j,LM)
           end if

           ! Compute condensation and thl, ql, qi
           call condensation_edmfA(upthv(iup,i,j), upqt(iup,i,j), ple(i,j,LM), &
                                   upthl(iup,i,j), upql(iup,i,j), upqi(iup,i,j), ice_ramp)
        end do
     else
        active_updraft(:,i,j) = .false.
     end if ! wthv > 0.
  end do
  end do

  ! Loop above surface
  do k = LM,2,-1
     km1 = k - 1
     kp1 = k + 1
     
     do j = 1,JM
     do i = 1,IM
        dzle = zle(i,j,km1) - zle(i,j,k)
        
        mfthvt = 0.
        mft    = 0.

        do iup = 1,numup
           if (active_updraft(iup,i,j)) then
              !
              ! Sample updraft statistics
              !

              ! Buoyancy flux at full levels for SHOC
              mfthvt_work = upa(iup,i,j)*upw(iup,i,j)*upthv(iup,i,j)
              mft_work    = upa(iup,i,j)*upw(iup,i,j)

              ! Statistics required for solver
              ltm = exfh(i,j,k)*( upthl(iup,i,j) - thli(i,j,k) )

              aw(i,j,k)    = aw(i,j,k)    + upa(iup,i,j)*upw(iup,i,j)
              aw2(i,j,k)   = aw2(i,j,k)   + upa(iup,i,j)*upw(iup,i,j)**2.
              ahl2(i,j,k)  = ahl2(i,j,k)  + upa(iup,i,j)*ltm**2.
              aqt2(i,j,k)  = aqt2(i,j,k)  + upa(iup,i,j)*( upqt(iup,i,j) - qti(i,j,k) )**2.
              aqthl(i,j,k) = aqthl(i,j,k) + upa(iup,i,j)*exfh(i,j,k)*( upthl(iup,i,j) - thli(i,j,k) )&
                                                                    *( upqt(iup,i,j) - qti(i,j,k) )
              aw3(i,j,k)    = aw3(i,j,k)   + upa(iup,i,j)*upw(iup,i,j)**3.
              aqt3(i,j,k)   = aqt3(i,j,k)  + upa(iup,i,j)*( upqt(iup,i,j) - qti(i,j,k) )**3.

              if ( edmf_implicit == 1 ) then
                 stmp =  mapl_cp*exfh(i,j,k)*upthl(iup,i,j) + mapl_grav*zle(i,j,k) + mapl_alhl*upql(iup,i,j) + mapl_alhs*upqi(iup,i,j) 
                        
                 awu(i,j,k)  = awu(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*upu(iup,i,j)
                 awv(i,j,k)  = awv(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*upv(iup,i,j)
                 aws(i,j,k)  = aws(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*stmp
                 awqv(i,j,k) = awqv(i,j,k) +  upa(iup,i,j)*upw(iup,i,j)&
                                             *( upqt(iup,i,j) - upql(iup,i,j) - upqi(iup,i,j) )
                 awql(i,j,k) = awql(i,j,k) + upa(iup,i,j)*upw(iup,i,j)*upql(iup,i,j)
                 awqi(i,j,k) = awqi(i,j,k) + upa(iup,i,j)*upw(iup,i,j)*upqi(iup,i,j)
              else
                 stmp =   mapl_cp*exfh(i,j,k)*( upthl(iup,i,j) - thli(i,j,k) ) &
                        + mapl_alhs*( upqi(iup,i,j) - qii(i,j,k) )             &
                        + mapl_alhl*( upql(iup,i,j) - qli(i,j,k) ) 

                 awu(i,j,k)  = awu(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*( upu(iup,i,j)  - ui(i,j,k) )
                 awv(i,j,k)  = awv(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*( upv(iup,i,j)  - vi(i,j,k) )
                 aws(i,j,k)  = aws(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*stmp
                 awqv(i,j,k) = awqv(i,j,k) +  upa(iup,i,j)*upw(iup,i,j)&
                                             *( upqt(iup,i,j)-upql(iup,i,j)-upqi(iup,i,j) - qvi(i,j,k))
                 awql(i,j,k) = awql(i,j,k) + upa(iup,i,j)*upw(iup,i,j)*( upql(iup,i,j) - qli(i,j,k) )
                 awqi(i,j,k) = awqi(i,j,k) + upa(iup,i,j)*upw(iup,i,j)*( upqi(iup,i,j) - qii(i,j,k) )
              end if

              ! Statistics required for MYNN-EDMF
              whl_mf(i,j,k)  = whl_mf(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*ltm
              wqt_mf(i,j,k)  = wqt_mf(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*( upqt(iup,i,j) - qti(i,j,k) )
              wthv_mf(i,j,k) = wthv_mf(i,j,k) + upa(iup,i,j)*upw(iup,i,j)*( upthv(iup,i,j) - thvi(i,j,k) )

              ! Dry-moist updraft statistics
              if ( upql(iup,i,j) > 0. .or. upqi(iup,i,j) > 0. ) then
                 edmfmoista(i,j,k)   = edmfmoista(i,j,k)   + upa(iup,i,j)
                 edmfmoistw(i,j,k)   = edmfmoistw(i,j,k)   + upa(iup,i,j)*upw(iup,i,j)
                 edmfmoistqt(i,j,k)  = edmfmoistqt(i,j,k)  + upa(iup,i,j)*upqt(iup,i,j)
                 edmfmoistthl(i,j,k) = edmfmoistthl(i,j,k) + upa(iup,i,j)*upthl(iup,i,j)
                 edmfmoistu(i,j,k)   = edmfmoistu(i,j,k)   + upa(iup,i,j)*upu(iup,i,j)
                 edmfmoistv(i,j,k)   = edmfmoistv(i,j,k)   + upa(iup,i,j)*upv(iup,i,j)
                 edmfmoistqc(i,j,k)  = edmfmoistqc(i,j,k)  + upa(iup,i,j)*( upql(iup,i,j) + upqi(iup,i,j) )
              else
                 edmfdrya(i,j,k)   = edmfdrya(i,j,k)   + upa(iup,i,j)
                 edmfdryw(i,j,k)   = edmfdryw(i,j,k)   + upa(iup,i,j)*upw(iup,i,j)
                 edmfdryqt(i,j,k)  = edmfdryqt(i,j,k)  + upa(iup,i,j)*upqt(iup,i,j)
                 edmfdrythl(i,j,k) = edmfdrythl(i,j,k) + upa(iup,i,j)*upthl(iup,i,j)
                 edmfdryu(i,j,k)   = edmfdryu(i,j,k)   + upa(iup,i,j)*upu(iup,i,j)
                 edmfdryv(i,j,k)   = edmfdryv(i,j,k)   + upa(iup,i,j)*upv(iup,i,j)
              end if

              !
              ! Compute fractional entrainment rate
              !

              ! Random number generator for stochastic entrainment
              call Poisson(1, 1, 1, 1, entf(iup,i,j,k), enti(iup,i,j,k), the_seed)

              if ( L0(i,j) > 0. ) then
                 ent(iup,i,j,k) = real(enti(iup,i,j,k))*ent0/dzle
!                 ent(iup,i,j,k) = entf(iup,i,j,k)*ent0/dzle ! test

                 ! increase entrainment if local minimum of thv
                 if ( ( thv(i,j,k) < thv(i,j,kp1) ) .and. ( thv(i,j,k) < thv(i,j,km1) ) ) then
                    ent(iup,i,j,k) = ent(iup,i,j,k) + 5.*ent0/L0(i,j)
                 end if
              else
                 ! negative L0 means 0 entrainment  
                 ent(:,i,j,:) = 0.
              end if

              ! sample entrainment
              E(i,j,k) = E(i,j,k) + rhoe(i,j,k)*upa(iup,i,j)*upw(iup,i,j)*ent(iup,i,j,k)

              !
              ! Integrate updraft equations
              !

              !
              EntExp  = exp(-ent(iup,i,j,k)*dzle)
              EntExpU = exp(-ent(iup,i,j,k)*dzle*EntWFac)

              ! Thermodynamic variables in updraft
              QTn  = qt(i,j,k)*( 1. - EntExp )  + upqt(iup,i,j)*EntExp
              THLn = thl(i,j,k)*( 1. - EntExp ) + upthl(iup,i,j)*EntExp
              Un   = u(i,j,k)*( 1. - EntExpU )  + upu(iup,i,j)*EntExpU
              Vn   = v(i,j,k)*( 1. - EntExpU )  + upv(iup,i,j)*EntExpU

              ! Condensation                 
              call condensation_edmf(QTn, THLn, ple(i,j,km1), THVn, QCn, wf, ice_ramp)

              ! Vertical velocity
              B  = mapl_grav*( 0.5*( THVn + upthv(iup,i,j) )/thv(i,j,k) - 1. )
              WP = Wb*ent(iup,i,j,k)
              if (WP == 0.) then
                 Wn2 = upw(iup,i,j)**2. + 2.*Wa*B*dzle
              else
                 EntW = exp(-2.*WP*dzle)
                 Wn2  = EntW*upw(iup,i,j)**2. + Wa*B/WP*( 1. - EntW )
              end if

              if ( Wn2 > 0. ) then
!                 UPW(K,I)   = sqrt(Wn2)
                 upw(iup,i,j)   = min( sqrt(Wn2), 5. ) ! npa  
                 upthv(iup,i,j) = THVn
                 upthl(iup,i,j) = THLn
                 upqt(iup,i,j)  = QTn
                 upql(iup,i,j)  = QCn*wf
                 upqi(iup,i,j)  = QCn*( 1. - wf )
                 upu(iup,i,j)   = Un
                 upv(iup,i,j)   = Vn

                 mfthvt = mfthvt + 0.5*( upa(iup,i,j)*upw(iup,i,j)*upthv(iup,i,j) + mfthvt_work )
                 mft    = mft    + 0.5*( upa(iup,i,j)*upw(iup,i,j) + mft_work )
              else
                 mfthvt = mfthvt + 0.5*mfthvt_work
                 mft    = mft    +  0.5*mft_work

                 active_updraft(iup,i,j) = .false.
              end if
           end if ! active_updraft(iup,i,j)
        end do ! iup = 1,numup

        buoyf(i,j,k) = buoyf(i,j,k) + exf(i,j,k)*( mfthvt - mft*thv(i,j,k) )

        ! Average updraft samples
        if ( edmfdrya(i,j,k) > 0. ) then
           edmfdryw(i,j,k)   = edmfdryw(i,j,k)/edmfdrya(i,j,k)
           edmfdryqt(i,j,k)  = edmfdryqt(i,j,k)/edmfdrya(i,j,k)
           edmfdrythl(i,j,k) = edmfdrythl(i,j,k)/edmfdrya(i,j,k)
           edmfdryu(i,j,k)   = edmfdryu(i,j,k)/edmfdrya(i,j,k)
           edmfdryv(i,j,k)   = edmfdryv(i,j,k)/edmfdrya(i,j,k)
        else
           edmfdryw(i,j,k)   = mapl_undef
           edmfdryqt(i,j,k)  = mapl_undef
           edmfdrythl(i,j,k) = mapl_undef
           edmfdryu(i,j,k)   = mapl_undef
           edmfdryv(i,j,k)   = mapl_undef
        end if

        if ( edmfmoista(i,j,k) > 0. ) then
           edmfmoistw(i,j,k)   = edmfmoistw(i,j,k)/edmfmoista(i,j,k)
           edmfmoistqt(i,j,k)  = edmfmoistqt(i,j,k)/edmfmoista(i,j,k)
           edmfmoistthl(i,j,k) = edmfmoistthl(i,j,k)/edmfmoista(i,j,k)
           edmfmoistu(i,j,k)   = edmfmoistu(i,j,k)/edmfmoista(i,j,k)
           edmfmoistv(i,j,k)   = edmfmoistv(i,j,k)/edmfmoista(i,j,k)
           edmfmoistqc(i,j,k)  = edmfmoistqc(i,j,k)/edmfmoista(i,j,k)
        else
           edmfmoistw(i,j,k)   = mapl_undef
           edmfmoistqt(i,j,k)  = mapl_undef
           edmfmoistthl(i,j,k) = mapl_undef
           edmfmoistu(i,j,k)   = mapl_undef
           edmfmoistv(i,j,k)   = mapl_undef
           edmfmoistqc(i,j,k)  = mapl_undef
        end if

        au(i,j,k) = edmfdrya(i,j,k) + edmfmoista(i,j,k)
        if ( au(i,j,k) > 0. ) then
           if ( edmfmoista(i,j,k) > 0. ) then
              wu(i,j,k)   = (  edmfdrya(i,j,k)*edmfdryw(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistw(i,j,k) )/au(i,j,k)
              thlu(i,j,k) = (  edmfdrya(i,j,k)*edmfdrythl(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistthl(i,j,k) )/au(i,j,k)
              qtu(i,j,k)  = (  edmfdrya(i,j,k)*edmfdryqt(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistqt(i,j,k) )/au(i,j,k)
           else
              wu(i,j,k)   = edmfdryw(i,j,k)
              thlu(i,j,k) = edmfdrythl(i,j,k)
              qtu(i,j,k)  = edmfdryqt(i,j,k)
           end if
        else
           wu(i,j,k)   = 0.
           thlu(i,j,k) = thli(i,j,k)
           qtu(i,j,k)  = qti(i,j,k)
        end if
        Mu(i,j,k) = rhoe(i,j,k)*au(i,j,k)*wu(i,j,k)

        ae(i,j,k) = ( 1. - edmfdrya(i,j,k) - edmfmoista(i,j,k) )*EDfac 

     end do ! i = 1,IM
     end do ! j = 1,JM
  end do ! k = LM,2,-1

  do k = 1,LM
     km1 = k - 1

     do j = 1,JM
     do i = 1,IM
        au_full   = 0.5*( au(i,j,k) + au(i,j,km1) )
        thlu_full = 0.5*( thlu(i,j,k) + thlu(i,j,km1) )
        qtu_full  = 0.5*( qtu(i,j,k) + qtu(i,j,km1) )

        if ( au_full > 0. ) then 
           hle(i,j,k) =   exf(i,j,k)*( thl(i,j,k) - au_full*thlu_full )/( 1. - au_full ) &
                        + (mapl_grav/mapl_cp)*z(i,j,k)
           qte(i,j,k) = ( qt(i,j,k) - au_full*qtu_full )/( 1. - au_full )
        else
           hle(i,j,k) = exf(i,j,k)*thl(i,j,k) + (mapl_grav/mapl_cp)*z(i,j,k)
           qte(i,j,k) = qt(i,j,k)
        end if

        D(i,j,k) = E(i,j,k) - ( Mu(i,j,km1) - Mu(i,j,k) )/( zle(i,j,km1) - zle(i,j,k) )

        buoyf(i,j,k) = buoyf(i,j,k) + exf(i,j,k)*( mfthvt - mft*thv(i,j,k) )
        
        mfw2(i,j,k)   = 0.5*( aw2(i,j,k)    + aw2(i,j,km1) )
        mfw3(i,j,k)   = 0.5*( aw3(i,j,k)    + aw3(i,j,km1) )
        mfhl2(i,j,k)  = 0.5*( ahl2(i,j,k)   + ahl2(i,j,km1) )
        mfqt2(i,j,k)  = 0.5*( aqt2(i,j,k)   + aqt2(i,j,km1) )
        mfqt3(i,j,k)  = 0.5*( aqt3(i,j,k)   + aqt3(i,j,km1) )
        mfwqt(i,j,k)  = 0.5*( wqt_mf(i,j,k) + wqt_mf(i,j,km1) )
        mfqthl(i,j,k) = 0.5*( aqthl(i,j,k)  + aqthl(i,j,km1) )
        mfwhl(i,j,k)  = 0.5*( whl_mf(i,j,k) + whl_mf(i,j,km1) )
     end do
     end do
  end do

end subroutine run_edmf


subroutine condensation_edmf(QT,THL,P,THV,QC,wf,ice_ramp)
!
! zero or one condensation for edmf: calculates THV and QC
!
use GEOS_UtilsMod, only : GEOS_Qsat

real,intent(in) :: QT,THL,P
real,intent(in) :: ice_ramp
real,intent(out):: THV,QC,wf


integer :: niter,i
real :: diff,exn,t,qs,qcold
 
! max number of iterations
niter=50
! minimum difference
diff=2.e-5

EXN=(P/mapl_p00)**mapl_kappa
QC=0. 
T=EXN*THL

do i=1,NITER
  T=EXN*THL+get_alhl(T,ice_ramp)/mapl_cp*QC
! qsat, p is in pascal
  QS=geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
  QCOLD=QC
  QC=max(0.5*QC+0.5*(QT-QS),0.)
if (abs(QC-QCOLD)<Diff) exit
enddo

T=EXN*THL+get_alhl(T,ice_ramp)/mapl_cp*QC
QS=geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
QC=max(QT-QS,0.)
THV=(THL+get_alhl(T,ice_ramp)/mapl_cp*QC/EXN)*(1.+(mapl_epsilon)*(QT-QC)-QC)
wf=water_f(T,ice_ramp)

end subroutine condensation_edmf


subroutine condensation_edmfA(THV,QT,P,THL,QL,QI,ice_ramp)
!
! zero or one condensation for edmf: calculates QL,QI from THV and QT
!

use GEOS_UtilsMod, only : GEOS_Qsat

real,intent(in) :: THV,QT,P
real,intent(in) :: ice_ramp
real,intent(out):: THL,QL,QI


integer :: niter,i
real :: diff,exn,t,qs,qcold,wf,qc
 
! max number of iterations
niter=50
! minimum difference
diff=2.e-5

EXN=(P/mapl_p00)**mapl_kappa
QC=0. 

T=EXN*THL

do i=1,NITER
   T=EXN*THV/(1.+mapl_epsilon*(QT-QC)-QC)
   QS=geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
   QCOLD=QC
   QC=max(0.5*QC+0.5*(QT-QS),0.)
if (abs(QC-QCOLD)<Diff) exit
enddo

 THL=(T-get_alhl(T,ice_ramp)/mapl_cp*QC)/EXN
 wf=water_f(T,ice_ramp)
 QL=QC*ice_ramp
 QI=QC*(1.-ice_ramp) 

end subroutine condensation_edmfA


function  get_alhl3(T,IM,JM,LM,iceramp)

real,dimension(IM,JM,LM) ::  T,get_alhl3
real :: iceramp
integer :: IM,JM,LM
integer :: I,J,L

do i=1,im
  do j=1,jm
    do l=1,lm
        get_alhl3(i,j,l)=get_alhl(T(i,j,l),iceramp)
    enddo
  enddo
enddo

end function get_alhl3

function get_alhl(T,iceramp)
   real :: T,get_alhl,iceramp,wf
    wf=water_f(T,iceramp)
    get_alhl=wf*mapl_alhl+(1.-wf)*mapl_alhs
end function get_alhl


function water_f(T,iceramp)
!
! computes water fraction
!
real ::T,iceramp,water_f,Tw
real :: Tmax,Tmin

  Tmax=0.
  Tmin=Tmax-abs(iceramp)
  Tw=T-273.16

! water fraction
  IF (Tw>Tmax) THEN
    water_f=1.
  ELSE IF (Tw<Tmin) THEN
    water_f=0.
  ELSE
    water_f=(Tw-Tmin)/(Tmax-Tmin);
  END IF

end function water_f


subroutine Poisson(istart,iend,jstart,jend,mu,POI,seed)

integer, intent(in) :: istart,iend,jstart,jend 
real,dimension(istart:iend,jstart:jend),intent(in) :: MU
integer, dimension(istart:iend,jstart:jend), intent(out) :: POI
integer :: sed_len
integer,dimension(2),  intent(in) :: seed

integer :: seed_len
integer :: i,j,idum,p
integer,allocatable :: theseed(:)

call random_seed(SIZE=seed_len)
allocate(theseed(seed_len))

theseed(1:2)=seed
! Gfortran uses longer seeds, so fill the rest with zero
if (seed_len > 2) theseed(3:) = seed(2)
 
 
call random_seed(put=theseed)


do i=istart,iend
 do j=jstart,jend
    poi(i,j)=poidev(mu(i,j),idum)
enddo
 enddo

end subroutine Poisson

subroutine Poisson_new(mu, poi, seed)

  real, intent(in)                  :: mu
  integer, dimension(2), intent(in) :: seed
  
  integer, intent(out) :: poi
  
  integer                            :: sed_len, seed_len, idum, p 
  integer, dimension(:), allocatable :: theseed

  call random_seed(SIZE=seed_len)
  allocate(theseed(seed_len))

  theseed(1:2) = seed
  ! Gfortran uses longer seeds, so fill the rest with zero
  if (seed_len > 2) theseed(3:) = seed(2)
  
  call random_seed(put=theseed)

  poi = poidev(mu,idum)

end subroutine Poisson_new



      FUNCTION poidev(xm,idum)
      INTEGER idum
      REAL poidev,xm,PI
      PARAMETER (PI=3.141592654)
!CU    USES gammln,ran1
      REAL alxm,em,g,oldm,sq,t,y
      SAVE alxm,g,oldm,sq
      DATA oldm /-1./
      if (xm.lt.12.)then
        if (xm.ne.oldm) then
          oldm=xm
          g=exp(-xm)
        endif
        em=-1
        t=1.
2       em=em+1.
        t=t*ran1(idum)
        if (t.gt.g) goto 2
      else
        if (xm.ne.oldm) then
          oldm=xm
          sq=sqrt(2.*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.)
        endif
1       y=tan(PI*ran1(idum))
        em=sq*y+xm
        if (em.lt.0.) goto 1
        em=int(em)
        t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
        if (ran1(idum).gt.t) goto 1
      endif
      poidev=em
      return
      END FUNCTION poidev

      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
     24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
     -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END FUNCTION gammln

      function ran1(idum)
      real ran1
      integer idum

      call random_number(ran1)
      end function ran1


end module edmf_mod

