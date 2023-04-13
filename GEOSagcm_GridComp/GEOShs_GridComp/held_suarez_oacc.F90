module hs_oacc_mod

  use MAPL_ConstantsMod, only: MAPL_CP, MAPL_GRAV, MAPL_KAPPA, MAPL_P00, MAPL_PI, MAPL_RGAS

contains

  subroutine held_suarez_oacc( &
       CPHI2, DISS, DTDT, DUDT, DVDT, &
       HFCN, P_I, PLE, SPHI2, TAUX, TAUY, T, &
       THEQ, T_EQ, U, V, &
       DAYLEN, DELH, DELV1, DT, FRICQ, FriendlyTemp, &
       FriendlyWind, GAM_D, GAM_I, IM, JM, LM, P_1, P_D, QMAX, &
       SIG1, TAUA, TAUF, TAUS, TSTRT, T0, compType, rank)

    use openacc

    implicit none

    integer, intent(in) :: FRICQ, IM, JM, LM, compType, rank
    logical, intent(in) :: FriendlyTemp, FriendlyWind

    real, dimension(:,:),   pointer,  intent(in)    :: CPHI2, HFCN, P_I, SPHI2
    real, dimension(:,:),   pointer,  intent(inout) :: DISS, TAUX, TAUY
    real, dimension(:,:,:), pointer,  intent(in)	:: PLE
    real, dimension(:,:,:), pointer,  intent(inout) :: DTDT, DUDT, DVDT, T, T_EQ, THEQ, U, V

    real, intent(in) :: DAYLEN, DELH, DELV1, DT, GAM_D, GAM_I, P_D, P_1
    real, intent(in) :: QMAX, SIG1, T0, TAUA, TAUF, TAUS, TSTRT

    real, dimension(IM,JM) :: DP, PL, UU, VV, VR, TE, DS, PII, F1, RR, DM, PK

    real :: DP_s, PL_s, UU_s, VV_s, VR_s, TE_s, DS_s, PII_s, F1_s, RR_s, DM_s, PK_s

    real, pointer, dimension(:,:) :: PS
    real, pointer, dimension(:,:) :: PT

    real :: t1, t2

    integer ::  L, I, J
    real    :: KA, KF, KS

    integer :: ngpus

    real :: CP, GRAV, KAPPA, P00, PI, RGAS

    CP = MAPL_CP
    GRAV = MAPL_GRAV
    KAPPA = MAPL_KAPPA
    P00 = MAPL_P00
    PI = MAPL_PI
    RGAS = MAPL_RGAS

    ! Initialize
    DP_s = 0.0
    PL_s = 0.0
    UU_s = 0.0
    VV_s = 0.0
    VR_s = 0.0
    TE_s = 0.0
    DS_s = 0.0
    PII_s = 0.0
    F1_s= 0.0
    RR_s = 0.0
    DM_s = 0.0
    PK_s = 0.0
    t1 = 0.0
    t2 = 0.0

    DP = 0.0
    PL = 0.0
    UU = 0.0
    VV = 0.0
    VR = 0.0
    TE = 0.0
    DS = 0.0
    PII = 0.0
    F1= 0.0
    RR = 0.0
    DM = 0.0
    PK = 0.0

    ngpus = acc_get_num_devices(acc_device_default)

    write(*,*) 'HS: rank =', rank, ', comp type: ', compType, ', Number of GPUS:', ngpus
    ! call acc_set_device_num(mod(rank,4),acc_device_nvidia)

    !write(*,*) 'From HS routine, rank =', rank, ": Number of GPUS = ", ngpus, ': This process is using', acc_get_device_num(acc_device_nvidia)

    ! If running with OpenMP, use these pointer assignments for PS and PT
    PS  => PLE(:,:,LM)
    PT  => PLE(:,:, 0)

    ! If running with OpenACC, use these pointer assignments for PS and PT
    !PS => PLE(:,:,LM+1)
    !PT => PLE(:,:,1)

    !if(associated(DISS)) DISS = 0.0
    !if(associated(TAUX)) TAUX = 0.0
    !if(associated(TAUY)) TAUY = 0.0

    KA   = 1.0/(DAYLEN*TAUA)
    KS   = 1.0/(DAYLEN*TAUS)
    KF   = 1.0/(DAYLEN*TAUF)

    PII =  P_D - (P_D - PT)*0.5*P_I

    !$acc data copyin(CPHI2, HFCN, P_I, SPHI2, PLE) &
    !$acc      copy(DISS, TAUX, TAUY) &
    !$acc      copy(DTDT, DUDT, DVDT, T, T_EQ, THEQ, U, V) &
    !$acc      copyin(PS, PT)

    ! write(*,*) 'Running HS Code with outer loop OpenMP/OpenACC Parallelism'
    call cpu_time(t1)

    !$acc parallel loop collapse(3) present(CPHI2, HFCN, P_I, SPHI2, PLE) &
    !$acc                           present(DISS, TAUX, TAUY) &
    !$acc                           present(DTDT, DUDT, DVDT, T, T_EQ, THEQ, U, V) &
    !$acc                           present(PS, PT) &
    !$acc                           private(DP_s, PL_s, UU_s, VV_s, VR_s, TE_s, DS_s, PII_s, F1_s, RR_s, DM_s, PK_s)
    do L = 1,LM
       do J = 1, JM
          do I = 1, IM

             ! If running with OpenMP, use these assignments for DP and PL
             DP_s  = (PLE(I,J,L)-PLE(I,J,L-1))
             PL_s  = (PLE(I,J,L)+PLE(I,J,L-1))*0.5

             ! If running with OpenACC, use these assignments for DP and PL
             !DP_s  = (PLE(I,J,L+1)-PLE(I,J,L))
             !PL_s  = (PLE(I,J,L+1)+PLE(I,J,L))*0.5
             DM_s  = DP_s / GRAV
             PK_s  = (PL_s/P00)**KAPPA

             ! H&S equilibrium temperature
             !----------------------------

             TE_s  = PK_s*( T0 - DELH*SPHI2(I,J) - DELV1*CPHI2(I,J)*log( PL_s/P00 ) )
             TE_s  = max( TE_s, TSTRT )

             ! Williamson Stratospheric modifications to equilibrium temperature
             ! -----------------------------------------------------------------

             if( PL_s < P_D ) then
                TE_s  = TSTRT*( min(1.0,PL_s/P_D)**(RGAS*GAM_D/GRAV)     &
                     + min(1.0,PL_s/PII_s)**(RGAS*GAM_I/GRAV) - 1 )
             endif

             !  Exports of equilibrium T and Theta
             !------------------------------------

             if(associated(T_EQ)) T_EQ(I,J,L) = TE_s
             ! T_EQ(I,J,L) = TE_s
             if(associated(THEQ)) THEQ(I,J,L) = TE_s/PK_s
             ! THEQ(I,J,L) = TE_s/PK_s

             ! Vertical structure of timescales in H&S.
             !---------------------------------------------

             F1_s  = max(0.0, ( (PL_s/PS(I,J))-SIG1 )/( 1.0-SIG1 ) )

             ! Atmospheric heating from H&S
             !-----------------------------

             RR_s = (KA + (KS-KA)*F1_s*CPHI2(I,J)**2) * (TE_s-T(I,J,L))

             if(associated(DTDT)) DTDT(I,J,L) = DP_s*RR_s
             ! DTDT(I,J,L) = DP_s*RR_s
             if(FriendlyTemp    ) T   (I,J,L) = T(I,J,L) + DT*RR_s

             ! Wind tendencies
             !----------------

             UU_s  = -U(I,J,L)*(F1_s*KF)
             VV_s  = -V(I,J,L)*(F1_s*KF)

             if(associated(DUDT)) DUDT(I,J,L) = UU_s
             ! DUDT(I,J,L) = UU_s
             if(associated(DVDT)) DVDT(I,J,L) = VV_s
             ! DVDT(I,J,L) = VV_s

             if(FriendlyWind) then
                U(I,J,L) = U(I,J,L) + DT*UU_s
                V(I,J,L) = V(I,J,L) + DT*VV_s
             end if

             !  Frictional heating from H&S drag
             !----------------------------------

             DS_s = U(I,J,L)*UU_s + V(I,J,L)*VV_s

             if(associated(DISS)) DISS(I,J) = DISS(I,J) - DS_s*DM_s
             ! DISS(I,J) = DISS(I,J) - DS_s*DM_s
             if(FRICQ /= 0) then
                if(associated(DTDT)) DTDT(I,J,L) = DTDT(I,J,L) - DS_s*(DP_s/CP  )
                ! DTDT(I,J,L) = DTDT(I,J,L) - DS_s*(DP_s/CP  )
                if(FriendlYTemp    ) T   (I,J,L) = T   (I,J,L) - DS_s*(DT/CP  )
             end if

             !  Surface stresses from vertically integrated H&S surface drag
             !--------------------------------------------------------------

             if(associated(TAUX)) TAUX(I,J) = TAUX(I,J) - UU_s*DM_s
             ! TAUX(I,J) = TAUX(I,J) - UU_s*DM_s
             if(associated(TAUY)) TAUY(I,J) = TAUY(I,J) - VV_s*DM_s
             ! TAUY(I,J) = TAUY(I,J) - VV_s*DM_s

             ! Localized heat source, if any
             !------------------------------

             if((associated(DTDT).or.FriendlyTemp) .and. QMAX/=0.0) then
                ! if(FriendlyTemp .and. QMAX/=0.0) then
                if(PL_s > P_1) then
                   VR_s = HFCN(I,J)*(QMAX/DAYLEN)*sin( PI*(PS(I,J)-PL_s)/(PS(I,J)-P_1) )
                else
                   VR_s = 0.
                endif

                if(associated(DTDT)) DTDT(I,J,L) = DTDT(I,J,L) + DP_s*VR_s
                ! DTDT(I,J,L) = DTDT(I,J,L) + DP_s*VR_s
                if(FriendlyTemp    ) T   (I,J,L) = T   (I,J,L) + DT*VR_s
             end if
          enddo
       enddo
    enddo
    !$acc end parallel loop
    call cpu_time(t2)
    !!$acc update host (DTDT, DUDT, DVDT, T_EQ, T, U, V)

    !$acc end data

    ! write(*,*) 'sum(T_EQ) = ', sum(T_EQ)

    write(*,*) 'Runtime = ', t2 - t1
    ! write(*,*) '***'

  end subroutine held_suarez_oacc

end module hs_oacc_mod
