module Process_Library_standalone

    use MAPL_ConstantsMod
    use GEOS_UtilsMod

    implicit none

    private

    ! In anvil/convective clouds
    real, parameter :: aT_ICE_ALL = 245.16
    real, parameter :: aT_ICE_MAX = 261.16
    real, parameter :: aICEFRPWR  = 2.0
    ! Over snow/ice SRF_TYPE = 2
    real, parameter :: iT_ICE_ALL = MAPL_TICE-40.0
    real, parameter :: iT_ICE_MAX = MAPL_TICE
    real, parameter :: iICEFRPWR  = 4.0
    ! Over Land     SRF_TYPE = 1
    real, parameter :: lT_ICE_ALL = 239.16
    real, parameter :: lT_ICE_MAX = 261.16
    real, parameter :: lICEFRPWR  = 2.0
    ! Over Oceans   SRF_TYPE = 0
    real, parameter :: oT_ICE_ALL = 238.16
    real, parameter :: oT_ICE_MAX = 263.16
    real, parameter :: oICEFRPWR  = 4.0

    ! parameters
    real, parameter :: EPSILON =  MAPL_H2OMW/MAPL_AIRMW
    real, parameter :: K_COND  =  2.4e-2    ! J m**-1 s**-1 K**-1
    real, parameter :: DIFFU   =  2.2e-5    ! m**2 s**-1
    ! LDRADIUS4
    real, parameter :: RHO_I   =  916.8     ! Density of ice crystal in kg/m^3
    real, parameter :: RHO_W   = 1000.0     ! Density of liquid water in kg/m^3
    real, parameter :: be      = 1./3. - 0.14
    real, parameter :: bx      = 100.* (3./(4.*MAPL_PI))**(1./3.) * 0.07*6.92
    ! combined constantc
    real, parameter :: cpbgrav = MAPL_CP/MAPL_GRAV
    real, parameter :: gravbcp = MAPL_GRAV/MAPL_CP
    real, parameter :: alhlbcp = MAPL_ALHL/MAPL_CP
    real, parameter :: alhfbcp = MAPL_ALHF/MAPL_CP
    real, parameter :: alhsbcp = MAPL_ALHS/MAPL_CP

    real, parameter :: R_AIR = 3.47e-3 !m3 Pa kg-1K-1

    real, parameter :: mapl_undef = 1.0e15  ! NOTE : This is the value pulled from MAPL_Mod

    public :: BUOYANCY2
!$acc declare create(gravbcp, alhlbcp,mapl_undef)
    contains

    subroutine BUOYANCY2( IM, JM, LM, T, Q, QS, DQS, DZ, ZLO, PLO, PS, SBCAPE, MLCAPE, MUCAPE, &
        SBCIN, MLCIN, MUCIN, BYNCY, LFC, LNB )

    ! Computes surface-based (SB), mixed-layer (ML) and most unstable (MU) versions 
    ! of CAPE and CIN. 

    integer,                intent(in)  :: IM, JM, LM
    real, dimension(:,:,:), intent(in)  :: T, Q, QS, DQS, DZ, ZLO, PLO
    real, dimension(:,:,:), intent(out) :: BYNCY
    real, pointer, dimension(:,:)       :: MLCAPE, MUCAPE, MLCIN, MUCIN
    real, dimension(:,:)                :: SBCAPE, SBCIN, LFC, LNB
    real, dimension(:,:),   intent(in)  :: PS

    real, dimension(IM,JM,LM)           :: Tve
    real, dimension(IM,JM)              :: MSEp, Qp, tmp1, tmp2
    integer, dimension(IM,JM)           :: Lev0

    integer :: I, J, L

    ! Tve   = T*(1.+MAPL_VIREPS*Q)
    ! BYNCY = MAPL_UNDEF
    ! MSEp  = 0.
    ! Qp    = 0.

!$acc data create(Tve, MSEp, Qp, tmp1, tmp2, Lev0)

    !$acc kernels
    Tve   = T*(1.+MAPL_VIREPS*Q)
    BYNCY = MAPL_UNDEF
    MSEp  = 0.
    Qp    = 0.
    !$acc end kernels

    ! Mixed-layer calculation. Parcel properties averaged over lowest 90 hPa
    if ( associated(MLCAPE) .and. associated(MLCIN) ) then
        !$acc kernels
        BYNCY = MAPL_UNDEF
        tmp1 = 0.
        Lev0 = LM
        do L = LM,1,-1
            where (PS-PLO(:,:,L).lt.90.) 
                MSEp = MSEp + (T(:,:,L) + gravbcp*ZLO(:,:,L) + alhlbcp*Q(:,:,L))*DZ(:,:,L) 
                Qp   = Qp   + Q(:,:,L)*DZ(:,:,L)
                tmp1 = tmp1 + DZ(:,:,L)
                Lev0 = L
            end where
            if (all(PS-PLO(:,:,L).gt.90.)) exit
        end do
        where (tmp1.gt.0.)   ! average
            MSEp = MSEp / tmp1
            Qp = Qp / tmp1
        end where
        !$acc end kernels
        call RETURN_CAPE_CIN_v2(ZLO, PLO, DZ,      & 
                                MSEp, Qp, Tve, QS, DQS,       &
                                MLCAPE, MLCIN, BYNCY, LFC, LNB, Lev0, PS, T, Q, MUCAPE, MUCIN, 0)
        !$acc kernels
       where (MLCAPE.le.0.)
            MLCAPE = MAPL_UNDEF
            MLCIN  = MAPL_UNDEF
       end where
       !$acc end kernels
    end if

    ! Most unstable calculation. Parcel in lowest 255 hPa with largest CAPE
    if ( associated(MUCAPE) .and. associated(MUCIN) ) then
        !$acc kernels
        MUCAPE = 0.
        MUCIN  = 0.
        BYNCY = MAPL_UNDEF
        LFC = MAPL_UNDEF
        LNB = MAPL_UNDEF
        !$acc end kernels
        ! do I = 1,IM
        !     do J = 1,JM
        !         do L = LM,1,-1

        do L = LM,1,-1
            do J = 1,JM  
                do I = 1,IM  
                    if (PS(I,J)-PLO(I,J,L).le.255.) then

                        MSEp(I,J) = T(I,J,L) + gravbcp*ZLO(I,J,L) + alhlbcp*Q(I,J,L)
                        Qp(I,J)   = Q(I,J,L)
                        call RETURN_CAPE_CIN( ZLO(I,J,1:L), PLO(I,J,1:L), DZ(I,J,1:L),      & 
                                            MSEp(I,J), Qp(I,J), Tve(I,J,1:L), QS(I,J,1:L), DQS(I,J,1:L),       &
                                            tmp1(I,J), tmp2(I,J), BYNCY(I,J,1:L), LFC(I,J), LNB(I,J) )
                        if (tmp1(I,J) .gt. MUCAPE(I,J)) then
                            MUCAPE(I,J) = tmp1(I,J)
                            MUCIN(I,J)  = tmp2(I,J)
                        end if
                    endif
                end do
            end do
        end do        

        where (MUCAPE.le.0.)
            MUCAPE = MAPL_UNDEF
            MUCIN  = MAPL_UNDEF
        end where
    end if

    ! Surface-based calculation
    MSEp = T(:,:,LM) + gravbcp*ZLO(:,:,LM) + alhlbcp*Q(:,:,LM)  ! parcel moist static energy
    Qp   = Q(:,:,LM)                                            ! parcel specific humidity
    do I = 1,IM
        do J = 1,JM
            call RETURN_CAPE_CIN( ZLO(I,J,:), PLO(I,J,:), DZ(I,J,:),      & 
                                    MSEp(I,J), Qp(I,J), Tve(I,J,:), QS(I,J,:), DQS(I,J,:),       &
                                    SBCAPE(I,J), SBCIN(I,J), BYNCY(I,J,:), LFC(I,J), LNB(I,J) )
        end do
    end do
    where (SBCAPE.le.0.)
        SBCAPE = MAPL_UNDEF
        SBCIN  = MAPL_UNDEF
    end where
!$acc end data
    end subroutine BUOYANCY2

    subroutine RETURN_CAPE_CIN( ZLO, PLO, DZ, MSEp, Qp, Tve, Qsate, DQS, CAPE, CIN, BYNCY, LFC, LNB )
        real,               intent(in)  :: MSEp, Qp
        real, dimension(:), intent(in)  :: ZLO, PLO, DZ, Tve, Qsate, DQS
        real,               intent(out) :: CAPE, CIN, LFC, LNB
        real, dimension(:), intent(out) :: BYNCY
    
        integer :: I, L, LM, KLNB, KLFC
        real    :: Qpnew, Tp, Tvp, Tlcl, Buoy, dq
        logical :: aboveLNB, aboveLFC, aboveLCL
    
        LM = size(ZLO,1)
    
        aboveLNB = .false.
        aboveLFC = .false.
    
        Qpnew = Qp
    
        CAPE = 0.
        CIN  = 0.
        BYNCY = 0.
        LFC = MAPL_UNDEF
        LNB = MAPL_UNDEF
    
        Tp = MSEp - gravbcp*ZLO(LM) - alhlbcp*Qp  ! initial parcel temp at source level LM
        Tlcl = find_tlcl( Tp, 100.*Qp/QSATE(LM) )
        aboveLCL = (Tp.lt.Tlcl)
      
        do L = LM-1,1,-1   ! start at level above source air
    
            ! determine parcel Qp, Tp      
            if ( .not. aboveLCL ) then
                Tp = Tp - gravbcp*(ZLO(L)-ZLO(L+1))                ! new parcel temperature w/o condensation 
                if (Tp.lt.Tlcl) then
                    Tp = Tp + gravbcp*(ZLO(L)-ZLO(L+1))             ! if cross LCL, revert Tp and go to aboveLCL below
                    aboveLCL = .true.
                end if
            end if
            if ( aboveLCL .and. Qpnew*alhlbcp.gt.0.01 ) then
                Tp = Tp - gravbcp*( ZLO(L)-ZLO(L+1) ) / ( 1.+alhlbcp*DQS(L) )     ! initial guess including condensation
                DO I = 1,10                                                       ! iterate until Qp=qsat(Tp)
                    dq = Qpnew - GEOS_QSAT( Tp, PLO(L) )
                    if (abs(dq*alhlbcp)<0.01) then
                        exit
                    end if
                    Tp = Tp + dq*alhlbcp/(1.+alhlbcp*DQS(L))
                    Qpnew = Qpnew - dq/(1.+alhlbcp*DQS(L))
                END DO
            end if
            Tp = MSEp - gravbcp*ZLO(L) - alhlbcp*Qpnew
            !  Qc = qp - qpnew.   ! condensate (not used for pseudoadiabatic ascent)
        
            Tvp = Tp*(1.+MAPL_VIREPS*Qpnew)              ! parcel virtual temp
            !  Tvp = Tp*(1.+0.61*Qpnew - Qc) ! condensate loading
        
            BYNCY(L) = MAPL_GRAV*(Tvp-Tve(L))/Tve(L)         ! parcel buoyancy
    
        end do
    
        ! if surface parcel immediately buoyant, scan upward to find first elevated
        ! B>0 level above a B<0 level, label it LFC.  If no such level, set LFC at surface.
        KLFC = LM
        KLNB = LM
        aboveLFC = .false.
        if (BYNCY(LM-1).gt.0.) then
            do L = LM-2,1,-1   ! scan up to find elevated LFC
                if (BYNCY(L).gt.0. .and. BYNCY(L+1).le.0.) then
                    KLFC = L
                    aboveLFC = .true.
                end if
                if (aboveLFC .and. BYNCY(L).lt.0. ) then 
                    KLNB = L
                    exit
                end if
            end do
        else   ! if surface parcel not immediately buoyant, LFC is first B>0 level
            do L = LM-1,1,-1
                if (BYNCY(L).gt.0. .and. .not.aboveLFC) then
                    KLFC = L
                    aboveLFC = .true.
                end if
                if (aboveLFC .and. BYNCY(L).lt.0.) then
                    KLNB = L
                    exit
                end if
            end do
        end if
        LFC = ZLO(KLFC)
        LNB = ZLO(KLNB)
    
        CIN = SUM( min(0.,BYNCY(KLFC:)*DZ(KLFC:)) )        ! define CIN as negative
    !    CAPE = SUM( max(0.,BYNCY(KLNB:KLFC)*DZ(KLNB:KLFC)) )
        CAPE = SUM( max(0.,BYNCY(:)*DZ(:)) )
    
    end subroutine RETURN_CAPE_CIN

    subroutine RETURN_CAPE_CIN_v2( ZLO, PLO, DZ, MSEp, Qp, Tve, Qsate, DQS, CAPE, CIN, BYNCY, LFC, LNB, &
                                   Lev0, PS, T, Q, MUCAPE, MUCIN, ctype )
        integer, dimension(:,:),   intent(in)  :: Lev0
        real,    dimension(:,:),   intent(in)  :: PS
        real,    dimension(:,:,:), intent(in)  :: ZLO, PLO, DZ, Tve, Qsate, DQS, T, Q
        real,    dimension(:,:),   intent(out) :: CAPE, CIN, LFC, LNB
        real,    dimension(:,:,:), intent(out) :: BYNCY
        real,    dimension(:,:),   intent(inout) :: MSEp, Qp, MUCAPE, MUCIN
        integer,                   intent(in)  :: ctype
    
        integer :: I, J, L, LL, IM, JM, LM, II, KLNB, KLFC
        real    :: Qpnew, Tp, Tvp, Tlcl, Buoy, dq
        logical :: aboveLNB, aboveLFC, aboveLCL
    
        IM = size(ZLO,1)
        JM = size(ZLO,2)

        if((ctype == 0) .or. (ctype == 2)) then
!$acc parallel loop gang vector collapse(2) &
!$acc          private(LM, aboveLNB, aboveLFC, Qpnew, &
!$acc                  Tp, Tlcl, aboveLCL, dq, Tvp, KLNB, KLFC)
            do J = 1,JM
                do I = 1, IM

                    if (ctype == 0) then
                        LM = Lev0(I,J)
                    else
                        LM = size(ZLO,3)
                    endif

                    aboveLNB = .false.
                    aboveLFC = .false.
                
                    Qpnew = Qp(I,J)
                    
                    CAPE(I,J) = 0.
                    CIN(I,J)  = 0.
                    BYNCY(I,J, 1:LM) = 0.
                    LFC(I,J) = MAPL_UNDEF
                    LNB(I,J) = MAPL_UNDEF
                
                    Tp = MSEp(I,J) - gravbcp*ZLO(I,J,LM) - alhlbcp*Qp(I,J)  ! initial parcel temp at source level LM
                    Tlcl = find_tlcl( Tp, 100.*Qp(I,J)/QSATE(I,J,LM) )
                    aboveLCL = (Tp.lt.Tlcl)
!$acc loop seq
                    do L = LM-1,1,-1   ! start at level above source air
                
                        ! determine parcel Qp, Tp      
                        if ( .not. aboveLCL ) then
                            Tp = Tp - gravbcp*(ZLO(I,J,L)-ZLO(I,J,L+1))                ! new parcel temperature w/o condensation 
                            if (Tp.lt.Tlcl) then
                                Tp = Tp + gravbcp*(ZLO(I,J,L)-ZLO(I,J,L+1))             ! if cross LCL, revert Tp and go to aboveLCL below
                                aboveLCL = .true.
                            end if
                        end if
                        if ( aboveLCL .and. Qpnew*alhlbcp.gt.0.01 ) then
                            Tp = Tp - gravbcp*( ZLO(I,J,L)-ZLO(I,J,L+1) ) / ( 1.+alhlbcp*DQS(I,J,L) )     ! initial guess including condensation
                            DO II = 1,10                                                       ! iterate until Qp=qsat(Tp)
                                dq = Qpnew - GEOS_QSAT( Tp, PLO(I,J,L) )
                                if (abs(dq*alhlbcp)<0.01) then
                                    exit
                                end if
                                Tp = Tp + dq*alhlbcp/(1.+alhlbcp*DQS(I,J,L))
                                Qpnew = Qpnew - dq/(1.+alhlbcp*DQS(I,J,L))
                            END DO
                        end if
                        Tp = MSEp(I,J) - gravbcp*ZLO(I,J,L) - alhlbcp*Qpnew
                        !  Qc = qp - qpnew.   ! condensate (not used for pseudoadiabatic ascent)
                    
                        Tvp = Tp*(1.+MAPL_VIREPS*Qpnew)              ! parcel virtual temp
                        !  Tvp = Tp*(1.+0.61*Qpnew - Qc) ! condensate loading
                    
                        BYNCY(I,J,L) = MAPL_GRAV*(Tvp-Tve(I,J,L))/Tve(I,J,L)         ! parcel buoyancy
                
                    end do
                
                    ! if surface parcel immediately buoyant, scan upward to find first elevated
                    ! B>0 level above a B<0 level, label it LFC.  If no such level, set LFC at surface.
                    KLFC = LM
                    KLNB = LM
                    aboveLFC = .false.
                    if (BYNCY(I,J,LM-1).gt.0.) then
                        do L = LM-2,1,-1   ! scan up to find elevated LFC
                            if (BYNCY(I,J,L).gt.0. .and. BYNCY(I,J,L+1).le.0.) then
                                KLFC = L
                                aboveLFC = .true.
                            end if
                            if (aboveLFC .and. BYNCY(I,J,L).lt.0. ) then 
                                KLNB = L
                                exit
                            end if
                        end do
                    else   ! if surface parcel not immediately buoyant, LFC is first B>0 level
                        do L = LM-1,1,-1
                            if (BYNCY(I,J,L).gt.0. .and. .not.aboveLFC) then
                                KLFC = L
                                aboveLFC = .true.
                            end if
                            if (aboveLFC .and. BYNCY(I,J,L).lt.0.) then
                                KLNB = L
                                exit
                            end if
                        end do
                    end if
                    LFC(I,J) = ZLO(I,J,KLFC)
                    LNB(I,J) = ZLO(I,J,KLNB)
                
                    CIN(I,J) = SUM( min(0.,BYNCY(I,J,KLFC:LM)*DZ(I,J,KLFC:LM)) )        ! define CIN as negative
                !    CAPE = SUM( max(0.,BYNCY(KLNB:KLFC)*DZ(KLNB:KLFC)) )
                    CAPE(I,J) = SUM( max(0.,BYNCY(I,J,1:LM)*DZ(I,J,1:LM)) )
                enddo
            enddo
!$acc end parallel
        else
            LM = size(Tve,3)
!$acc parallel loop gang vector collapse(2) &
!$acc          private(aboveLNB, aboveLFC, Qpnew, &
!$acc                  Tp, Tlcl, aboveLCL, dq, Tvp, KLNB, KLFC, LL)
            do J = 1,IM
                do I = 1, IM
!$acc loop seq
                    do LL = LM, 1, -1
                        if(PS(I,J) - PLO(I,J,LL) .gt. 255.) exit
                        MSEp(I,J) = T(I,J,LL) + gravbcp*ZLO(I,J,LL) + alhlbcp*Q(I,J,LL)
                        Qp(I,J)   = Q(I,J,LL)
                        aboveLNB = .false.
                        aboveLFC = .false.
                    
                        Qpnew = Qp(I,J)
                        
                        CAPE(I,J) = 0.
                        CIN(I,J)  = 0.
                        BYNCY(I,J, 1:LL) = 0.
                        LFC(I,J) = MAPL_UNDEF
                        LNB(I,J) = MAPL_UNDEF
                    
                        Tp = MSEp(I,J) - gravbcp*ZLO(I,J,LL) - alhlbcp*Qp(I,J)  ! initial parcel temp at source level LM
                        Tlcl = find_tlcl( Tp, 100.*Qp(I,J)/QSATE(I,J,LL) )
                        aboveLCL = (Tp.lt.Tlcl)
                    
                        do L = LL-1,1,-1   ! start at level above source air
                    
                            ! determine parcel Qp, Tp      
                            if ( .not. aboveLCL ) then
                                Tp = Tp - gravbcp*(ZLO(I,J,L)-ZLO(I,J,L+1))                ! new parcel temperature w/o condensation 
                                if (Tp.lt.Tlcl) then
                                    Tp = Tp + gravbcp*(ZLO(I,J,L)-ZLO(I,J,L+1))             ! if cross LCL, revert Tp and go to aboveLCL below
                                    aboveLCL = .true.
                                end if
                            end if
                            if ( aboveLCL .and. Qpnew*alhlbcp.gt.0.01 ) then
                                Tp = Tp - gravbcp*( ZLO(I,J,L)-ZLO(I,J,L+1) ) / ( 1.+alhlbcp*DQS(I,J,L) )     ! initial guess including condensation
                                DO II = 1,10                                                       ! iterate until Qp=qsat(Tp)
                                    dq = Qpnew - GEOS_QSAT( Tp, PLO(I,J,L) )
                                    if (abs(dq*alhlbcp)<0.01) then
                                        exit
                                    end if
                                    Tp = Tp + dq*alhlbcp/(1.+alhlbcp*DQS(I,J,L))
                                    Qpnew = Qpnew - dq/(1.+alhlbcp*DQS(I,J,L))
                                END DO
                            end if
                            Tp = MSEp(I,J) - gravbcp*ZLO(I,J,L) - alhlbcp*Qpnew
                            !  Qc = qp - qpnew.   ! condensate (not used for pseudoadiabatic ascent)
                        
                            Tvp = Tp*(1.+MAPL_VIREPS*Qpnew)              ! parcel virtual temp
                            !  Tvp = Tp*(1.+0.61*Qpnew - Qc) ! condensate loading
                        
                            BYNCY(I,J,L) = MAPL_GRAV*(Tvp-Tve(I,J,L))/Tve(I,J,L)         ! parcel buoyancy
                    
                        end do
                    
                        ! if surface parcel immediately buoyant, scan upward to find first elevated
                        ! B>0 level above a B<0 level, label it LFC.  If no such level, set LFC at surface.
                        KLFC = LL
                        KLNB = LL
                        aboveLFC = .false.
                        if (BYNCY(I,J,LL-1).gt.0.) then
                            do L = LL-2,1,-1   ! scan up to find elevated LFC
                                if (BYNCY(I,J,L).gt.0. .and. BYNCY(I,J,L+1).le.0.) then
                                    KLFC = L
                                    aboveLFC = .true.
                                end if
                                if (aboveLFC .and. BYNCY(I,J,L).lt.0. ) then 
                                    KLNB = L
                                    exit
                                end if
                            end do
                        else   ! if surface parcel not immediately buoyant, LFC is first B>0 level
                            do L = LL-1,1,-1
                                if (BYNCY(I,J,L).gt.0. .and. .not.aboveLFC) then
                                    KLFC = L
                                    aboveLFC = .true.
                                end if
                                if (aboveLFC .and. BYNCY(I,J,L).lt.0.) then
                                    KLNB = L
                                    exit
                                end if
                            end do
                        end if
                        LFC(I,J) = ZLO(I,J,KLFC)
                        LNB(I,J) = ZLO(I,J,KLNB)
                    
                        CIN(I,J) = SUM( min(0.,BYNCY(I,J,KLFC:LL)*DZ(I,J,KLFC:LL)) )        ! define CIN as negative
                    !    CAPE = SUM( max(0.,BYNCY(KLNB:KLFC)*DZ(KLNB:KLFC)) )
                        CAPE(I,J) = SUM( max(0.,BYNCY(I,J,1:LL)*DZ(I,J,1:LL)) )
                    enddo
                    if (CAPE(I,J) .gt. MUCAPE(I,J)) then
                        MUCAPE(I,J) = CAPE(I,J)
                        MUCIN(I,J)  = CIN(I,J)
                    end if
                enddo
            enddo
!$acc end parallel
        endif
    
    end subroutine RETURN_CAPE_CIN_v2

    FUNCTION FIND_TLCL ( tk, rh ) result( tlcl )
!$acc routine seq
        ! Description:                                                            
        !    This function calculates the temperature of a parcel of air would have
        !    if lifed dry adiabatically to it's lifting condensation level (lcl).  
        ! References:                                                              
        !    Bolton (1980), Monthly Weather Review, pg. 1048, Eq. 22
            IMPLICIT NONE
            REAL, INTENT ( IN ) :: tK   !~ Temperature ( K )
            REAL, INTENT ( IN ) :: rh   !~ Relative Humidity ( % )
            REAL                :: tlcl
            REAL :: denom, term1, term2
            term1 = 1.0 / ( tK - 55.0 )
            term2 = ( LOG (max(0.1,rh)/100.0)  / 2840.0 )
            denom = term1 - term2
            tlcl = ( 1.0 / denom ) + 55.0 
    END FUNCTION FIND_TLCL
end module

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”
  
! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”
  
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at
  
! http://www.apache.org/licenses/LICENSE-2.0
  
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.