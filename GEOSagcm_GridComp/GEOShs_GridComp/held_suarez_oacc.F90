module hs_oacc_mod

  use MAPL_ConstantsMod, only: MAPL_CP, MAPL_GRAV, MAPL_KAPPA, MAPL_P00, MAPL_PI, MAPL_RGAS
  use iso_c_binding, only: c_int, c_float, c_loc, c_f_pointer

contains

  subroutine hs_buffer_layer( &
       CPHI2, DISS, DTDT, DUDT, DVDT, &
       HFCN, P_I, PLE, SPHI2, TAUX, TAUY, T, &
       THEQ, T_EQ, U, V, &
       DAYLEN, DELH, DELV1, DT, FRICQ, FriendlyTemp, &
       FriendlyWind, GAM_D, GAM_I, IM, JM, LM, P_1, P_D, QMAX, &
       SIG1, TAUA, TAUF, TAUS, TSTRT, T0, compType, rank)

    implicit none

    integer, intent(in) :: FRICQ, IM, JM, LM, compType, rank
    logical, intent(in) :: FriendlyTemp, FriendlyWind

    real, dimension(IM,JM), target,		   intent(in)    :: CPHI2, HFCN, P_I, SPHI2
    real, dimension(IM,JM), target,		   intent(inout) :: DISS, TAUX, TAUY
    real, dimension(IM, JM, 0:LM), target, intent(in)    :: PLE
    real, dimension(IM, JM, LM),  target,  intent(inout) :: DTDT, DUDT, DVDT, T, T_EQ, THEQ, U, V

    real, intent(in) :: DAYLEN, DELH, DELV1, DT, GAM_D, GAM_I, P_D, P_1
    real, intent(in) :: QMAX, SIG1, T0, TAUA, TAUF, TAUS, TSTRT

    real, pointer, dimension(:,:) :: DISS_P, TAUX_P, TAUY_P
    real, pointer, dimension(:,:,:) :: T_EQ_P, THEQ_P, DTDT_P, DUDT_P, DVDT_P, PLE_P, T_P, U_P, V_P
    real, pointer, dimension(:,:) :: CPHI2_P, HFCN_P, P_I_P, SPHI2_P

    CPHI2_P => CPHI2
    HFCN_P => HFCN
    P_I_P => P_I
    SPHI2_P => SPHI2
    DISS_P => DISS
    TAUX_P => TAUX
    TAUY_P => TAUY
    PLE_P => PLE
    DTDT_P => DTDT
    DUDT_P => DUDT
    DVDT_P => DVDT
    T_P => T
    T_EQ_P => T_EQ
    THEQ_P => THEQ
    U_P => U
    V_P => V

    ! write(*,*) 'Going into hs_oacc'

    call held_suarez_oacc( &
         CPHI2_P, DISS_P, DTDT_P, DUDT_P, DVDT_P, &
         HFCN_P, P_I_P, PLE_P, SPHI2_P, TAUX_P, TAUY_P, T_P, &
         THEQ_P, T_EQ_P, U_P, V_P, &
         DAYLEN, DELH, DELV1, DT, FRICQ, FriendlyTemp, &
         FriendlyWind, GAM_D, GAM_I, IM, JM, LM, P_1, P_D, QMAX, &
         SIG1, TAUA, TAUF, TAUS, TSTRT, T0, compType, rank)

  end subroutine hs_buffer_layer

  subroutine c_held_suarez_oacc( &
       CPHI2, DISS, DTDT, DUDT, DVDT, &
       HFCN, P_I, PLE, SPHI2, TAUX, TAUY, T, &
       THEQ, T_EQ, U, V, &
       DAYLEN, DELH, DELV1, DT, FRICQ, FriendlyTemp, &
       FriendlyWind, GAM_D, GAM_I, IM, JM, LM, P_1, P_D, QMAX, &
       SIG1, TAUA, TAUF, TAUS, TSTRT, T0, compType, rank) bind(C)

    implicit none

    integer(C_INT), value :: FRICQ, IM, JM, LM, compType, rank
    logical(C_INT), value :: FriendlyTemp, FriendlyWind

    real(C_FLOAT), dimension(IM*JM), target :: CPHI2, HFCN, P_I, SPHI2
    real(C_FLOAT), dimension(IM*JM), target :: DISS, TAUX, TAUY
    real(C_FLOAT), dimension(IM*JM*(LM+1)), target :: PLE
    real(C_FLOAT), dimension(IM*JM*LM), target :: DTDT, DUDT, DVDT, T, T_EQ, THEQ, U, V

    real, dimension(:,:), pointer :: CPHI2_p, HFCN_p, P_I_p, SPHI2_p
    real, dimension(:,:), pointer :: DISS_p, TAUX_p, TAUY_p
    real, dimension(:,:,:), pointer :: PLE_p
    real, dimension(:,:,:), pointer :: DTDT_p, DUDT_p, DVDT_p, T_p, T_EQ_p, THEQ_p, U_p, V_p

    real(C_FLOAT), value :: DAYLEN, DELH, DELV1, DT, GAM_D, GAM_I, P_D, P_1
    real(C_FLOAT), value :: QMAX, SIG1, T0, TAUA, TAUF, TAUS, TSTRT

    call c_f_pointer(C_LOC(CPHI2), CPHI2_p, [IM, JM])
    call c_f_pointer(C_LOC(HFCN),  HFCN_p,  [IM, JM])
    call c_f_pointer(C_LOC(P_I),   P_I_p,   [IM, JM])
    call c_f_pointer(C_LOC(SPHI2), SPHI2_p, [IM, JM])
    call c_f_pointer(C_LOC(DISS),  DISS_p,  [IM, JM])
    call c_f_pointer(C_LOC(TAUX),  TAUX_p,  [IM, JM])
    call c_f_pointer(C_LOC(TAUY),  TAUY_p,  [IM, JM])
    ! call c_f_pointer(C_LOC(PLE),   PLE_p,   [IM, JM, 0:LM])
    call c_f_pointer(C_LOC(PLE),   PLE_p,   [IM, JM, LM+1])
    call c_f_pointer(C_LOC(DTDT),  DTDT_p,  [IM, JM, LM])
    call c_f_pointer(C_LOC(DUDT),  DUDT_p,  [IM, JM, LM])
    call c_f_pointer(C_LOC(DVDT),  DVDT_p,  [IM, JM, LM])
    call c_f_pointer(C_LOC(T),     T_p,     [IM, JM, LM])
    call c_f_pointer(C_LOC(T_EQ),  T_EQ_p,  [IM, JM, LM])
    call c_f_pointer(C_LOC(THEQ),  THEQ_p,  [IM, JM, LM])
    call c_f_pointer(C_LOC(U),     U_p,     [IM, JM, LM])
    call c_f_pointer(C_LOC(V),     V_p,     [IM, JM, LM])

    call held_suarez_oacc( &
         CPHI2_P, DISS_P, DTDT_P, DUDT_P, DVDT_P, &
         HFCN_P, P_I_P, PLE_P, SPHI2_P, TAUX_P, TAUY_P, T_P, &
         THEQ_P, T_EQ_P, U_P, V_P, &
         DAYLEN, DELH, DELV1, DT, FRICQ, FriendlyTemp, &
         FriendlyWind, GAM_D, GAM_I, IM, JM, LM, P_1, P_D, QMAX, &
         SIG1, TAUA, TAUF, TAUS, TSTRT, T0, compType, rank)

  end subroutine c_held_suarez_oacc

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

    ngpus = acc_get_num_devices(acc_device_nvidia)

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

    !$acc data copyin(CPHI2, HFCN, P_I, SPHI2) &
    !$acc      copyin(DISS, TAUX, TAUY) &
    !$acc      copyin(PLE,DTDT, DUDT, DVDT, T, T_EQ, THEQ, U, V) &
    !$acc      copyin(DAYLEN, DELH, DELV1, DT, GAM_D, GAM_I, P_D, P_1) &
    !$acc      copyin(QMAX, SIG1, T0, TAUA, TAUF, TAUS, TSTRT) &
    !$acc      copyin(DP, PL, UU, VV, VR, TE, DS, PII, F1, RR, DM, PK) &
    !$acc      copyin(DP_s, PL_s, UU_s, VV_s, VR_s, TE_s, DS_s, PII_s, F1_s, RR_s, DM_s, PK_s) &
    !$acc      copyin(PS, PT) &
    !$acc      copyin(KA, KF, KS) &
    !$acc      copyin(CP, GRAV, KAPPA, P00, PI, RGAS) &
    !$acc      copyin(FRICQ, IM, JM, LM, compType, FriendlyTemp, FriendlyWind)

    if (compType == 0) then
       ! write(*,*) 'Running original HS Code'
       call cpu_time(t1)

       LEVELS: do L = 1,LM

          ! If running with OpenMP, use these assignments for DP and PL
          DP  = (PLE(:,:,L)-PLE(:,:,L-1))
          PL  = (PLE(:,:,L)+PLE(:,:,L-1))*0.5

          ! If running with OpenACC, use these assignments for DP and PL
          !DP  = (PLE(:,:,L+1)-PLE(:,:,L))
          !PL  = (PLE(:,:,L+1)+PLE(:,:,L))*0.5
          DM  = DP / GRAV
          PK  = (PL/P00)**KAPPA 

          ! H&S equilibrium temperature
          !----------------------------

          TE  = PK*( T0 - DELH*SPHI2 - DELV1*CPHI2*log( PL/P00 ) )
          TE  = max( TE, TSTRT )

          ! Williamson Stratospheric modifications to equilibrium temperature
          ! -----------------------------------------------------------------

          where( PL < P_D )
             TE  = TSTRT*( min(1.0,PL/P_D)**(RGAS*GAM_D/GRAV)     &
                  + min(1.0,PL/PII)**(RGAS*GAM_I/GRAV) - 1 )
          end where

          !  Exports of equilibrium T and Theta
          !------------------------------------

          if(associated(T_EQ)) T_EQ(:,:,L) = TE
          ! T_EQ(:,:,L) = TE
!!!!if(associated(THEQ)) THEQ(:,:,L) = TE/PK
          ! THEQ(:,:,L) = TE/PK

          ! Vertical structure of timescales in H&S.
          !---------------------------------------------

          F1  = max(0.0, ( (PL/PS)-SIG1 )/( 1.0-SIG1 ) )

          ! Atmospheric heating from H&S
          !-----------------------------

          RR = (KA + (KS-KA)*F1*CPHI2**2) * (TE-T(:,:,L))

          if(associated(DTDT)) DTDT(:,:,L) =            DP*RR
          ! DTDT(:,:,L) = DP*RR
          if(FriendlyTemp    ) T   (:,:,L) = T(:,:,L) + DT*RR

          ! Wind tendencies
          !----------------

          UU  = -U(:,:,L)*(F1*KF)
          VV  = -V(:,:,L)*(F1*KF)

          if(associated(DUDT)) DUDT(:,:,L) = UU
          ! DUDT(:,:,L) = UU
          if(associated(DVDT)) DVDT(:,:,L) = VV
          ! DVDT(:,:,L) = VV

          if(FriendlyWind) then
             U(:,:,L) = U(:,:,L) + DT*UU
             V(:,:,L) = V(:,:,L) + DT*VV
          end if

          !  Frictional heating from H&S drag
          !----------------------------------

          DS = U(:,:,L)*UU + V(:,:,L)*VV

!!!!if(associated(DISS)) DISS           = DISS        - DS*DM
          ! DISS = DISS - DS*DM
          if(FRICQ /= 0) then
             if(associated(DTDT)) DTDT(:,:,L) = DTDT(:,:,L) - DS*(DP/CP  )
             ! DTDT(:,:,L) = DTDT(:,:,L) - DS*(DP/CP  )
             if(FriendlYTemp    ) T   (:,:,L) = T   (:,:,L) - DS*(DT/CP  )
          end if

          !  Surface stresses from vertically integrated H&S surface drag
          !--------------------------------------------------------------

!!!if(associated(TAUX)) TAUX = TAUX - UU*DM
          ! TAUX = TAUX - UU*DM
!!!if(associated(TAUY)) TAUY = TAUY - VV*DM
          ! TAUY = TAUY - VV*DM

          ! Localized heat source, if any
          !------------------------------

          if((associated(DTDT).or.FriendlyTemp) .and. QMAX/=0.0) then
             ! if(FriendlyTemp .and. QMAX/=0.0) then
             where(PL > P_1)
                VR = HFCN*(QMAX/DAYLEN)*sin( PI*(PS-PL)/(PS-P_1) )
             elsewhere
                VR = 0.
             end where

             if(associated(DTDT)) DTDT(:,:,L) = DTDT(:,:,L) + DP*VR 
             ! DTDT(:,:,L) = DTDT(:,:,L) + DP*VR 
             if(FriendlyTemp    ) T   (:,:,L) = T   (:,:,L) + DT*VR
          end if

       enddo LEVELS
       call cpu_time(t2)

    else if (compType == 1) then
       ! write(*,*) 'Running HS Code with inner loop OpenMP/OpenACC Parallelism'
       call cpu_time(t1)

       do L = 1,LM		
          !$acc parallel loop collapse(2) present(CPHI2, HFCN, P_I, SPHI2) &
          !$acc                           present(DISS, TAUX, TAUY) &
          !$acc                           present(PLE,DTDT, DUDT, DVDT, T, T_EQ, THEQ, U, V) &
          !$acc                           present(DAYLEN, DELH, DELV1, DT, GAM_D, GAM_I, P_D, P_1) &
          !$acc                           present(QMAX, SIG1, T0, TAUA, TAUF, TAUS, TSTRT) &
          !$acc                           present(DP, PL, UU, VV, VR, TE, DS, PII, F1, RR, DM, PK) &
          !$acc                           present(PS, PT, KA, KF, KS) &
          !$acc                           present(CP, GRAV, KAPPA, P00, PI, RGAS) &
          !$acc                           present(FRICQ, IM, JM, LM, compType, FriendlyTemp, FriendlyWind)
          !$omp parallel do collapse(2) default(shared)
          do J = 1, JM
             do I = 1, IM

                ! If running with OpenMP, use these assignments for DP and PL
                DP(I,J)  = (PLE(I,J,L)-PLE(I,J,L-1))
                PL(I,J)  = (PLE(I,J,L)+PLE(I,J,L-1))*0.5

                ! If running with OpenACC, use these assignments for DP and PL
                !DP(I,J)  = (PLE(I,J,L+1)-PLE(I,J,L))
                !PL(I,J)  = (PLE(I,J,L+1)+PLE(I,J,L))*0.5

                DM(I,J)  = DP(I,J) / GRAV
                PK(I,J)  = (PL(I,J)/P00)**KAPPA 

                ! H&S equilibrium temperature
                !----------------------------

                TE(I,J)  = PK(I,J)*( T0 - DELH*SPHI2(I,J) - DELV1*CPHI2(I,J)*log( PL(I,J)/P00 ) )
                TE(I,J)  = max( TE(I,J), TSTRT )

                ! Williamson Stratospheric modifications to equilibrium temperature
                ! -----------------------------------------------------------------

                if( PL(I,J) < P_D ) then
                   TE(I,J)  = TSTRT*( min(1.0,PL(I,J)/P_D)**(RGAS*GAM_D/GRAV)     &
                        + min(1.0,PL(I,J)/PII(I,J))**(RGAS*GAM_I/GRAV) - 1 )
                endif

                !  Exports of equilibrium T and Theta
                !------------------------------------

                if(associated(T_EQ)) T_EQ(I,J,L) = TE(I,J)
                ! T_EQ(I,J,L) = TE(I,J)
                if(associated(THEQ)) THEQ(I,J,L) = TE(I,J)/PK(I,J)
                ! THEQ(I,J,L) = TE(I,J)/PK(I,J)

                ! Vertical structure of timescales in H&S.
                !---------------------------------------------

                F1(I,J)  = max(0.0, ( (PL(I,J)/PS(I,J))-SIG1 )/( 1.0-SIG1 ) )

                ! Atmospheric heating from H&S
                !-----------------------------

                RR(I,J) = (KA + (KS-KA)*F1(I,J)*CPHI2(I,J)**2) * (TE(I,J)-T(I,J,L))

                if(associated(DTDT)) DTDT(I,J,L) = DP(I,J)*RR(I,J)
                ! DTDT(I,J,L) = DP(I,J)*RR(I,J)
                if(FriendlyTemp    ) T   (I,J,L) = T(I,J,L) + DT*RR(I,J)

                ! Wind tendencies
                !----------------

                UU(I,J)  = -U(I,J,L)*(F1(I,J)*KF)
                VV(I,J)  = -V(I,J,L)*(F1(I,J)*KF)

                if(associated(DUDT)) DUDT(I,J,L) = UU(I,J)
                ! DUDT(I,J,L) = UU(I,J)
                if(associated(DVDT)) DVDT(I,J,L) = VV(I,J)
                ! DVDT(I,J,L) = VV(I,J)

                if(FriendlyWind) then
                   U(I,J,L) = U(I,J,L) + DT*UU(I,J)
                   V(I,J,L) = V(I,J,L) + DT*VV(I,J)
                end if

                !  Frictional heating from H&S drag
                !----------------------------------

                DS(I,J) = U(I,J,L)*UU(I,J) + V(I,J,L)*VV(I,J)

                if(associated(DISS)) DISS(I,J) = DISS(I,J) - DS(I,J)*DM(I,J)
                ! DISS(I,J) = DISS(I,J) - DS(I,J)*DM(I,J)
                if(FRICQ /= 0) then
                   if(associated(DTDT)) DTDT(I,J,L) = DTDT(I,J,L) - DS(I,J)*(DP(I,J)/CP  )
                   ! DTDT(I,J,L) = DTDT(I,J,L) - DS(I,J)*(DP(I,J)/CP  )
                   if(FriendlYTemp    ) T   (I,J,L) = T   (I,J,L) - DS(I,J)*(DT/CP  )
                end if

                !  Surface stresses from vertically integrated H&S surface drag
                !--------------------------------------------------------------

                if(associated(TAUX)) TAUX(I,J) = TAUX(I,J) - UU(I,J)*DM(I,J)
                ! TAUX(I,J) = TAUX(I,J) - UU(I,J)*DM(I,J)
                if(associated(TAUY)) TAUY(I,J) = TAUY(I,J) - VV(I,J)*DM(I,J)
                ! TAUY(I,J) = TAUY(I,J) - VV(I,J)*DM(I,J)

                ! Localized heat source, if any
                !------------------------------

                if((associated(DTDT).or.FriendlyTemp) .and. QMAX/=0.0) then
                   ! if(FriendlyTemp .and. QMAX/=0.0) then
                   if(PL(I,J) > P_1) then
                      VR(I,J) = HFCN(I,J)*(QMAX/DAYLEN)*sin( PI*(PS(I,J)-PL(I,J))/(PS(I,J)-P_1) )
                   else
                      VR(I,J) = 0.
                   endif

                   if(associated(DTDT)) DTDT(I,J,L) = DTDT(I,J,L) + DP(I,J)*VR(I,J)
                   ! DTDT(I,J,L) = DTDT(I,J,L) + DP(I,J)*VR(I,J)
                   if(FriendlyTemp    ) T   (I,J,L) = T   (I,J,L) + DT*VR(I,J)
                end if
             enddo
          enddo
          !$omp end parallel do
          !$acc end parallel loop
       enddo
       call cpu_time(t2)
       !$acc update host (DTDT, DUDT, DVDT, T_EQ, T, U, V)
    else
       ! write(*,*) 'Running HS Code with outer loop OpenMP/OpenACC Parallelism'
       call cpu_time(t1)

       !$omp parallel do collapse(3) default(shared) &
       !$omp private (DP_s, PL_s, UU_s, VV_s, VR_s, TE_s, DS_s, PII_s, F1_s, RR_s, DM_s, PK_s)
       !$acc parallel loop collapse(3) present(DP, PLE, PL, DM, GRAV, PK, KAPPA) &
       !$acc                           present(TE, T0, DELH, SPHI2, DELV1, CPHI2, P00) &
       !$acc                           present(TSTRT, P_D, RGAS, GAM_D) &
       !$acc                           present(PII, GAM_I, T_EQ, THEQ, F1, PS, SIG1, RR) &
       !$acc                           present(FriendlyTemp, T, DT, UU, VV, U, V, DS, DISS) &
       !$acc                           present(FRICQ, DUDT, DVDT, DTDT, CP, TAUX, TAUY, P_1, VR, HFCN) &
       !$acc                           present(QMAX, DAYLEN, PI, KF, KA, KS) &
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
       !$omp end parallel do
       !$acc end parallel loop
       call cpu_time(t2)
       !$acc update host (DTDT, DUDT, DVDT, T_EQ, T, U, V)

    endif
    !$acc end data

    ! write(*,*) 'sum(T_EQ) = ', sum(T_EQ)

    write(*,*) 'Runtime = ', t2 - t1
    ! write(*,*) '***'

  end subroutine held_suarez_oacc

end module hs_oacc_mod
