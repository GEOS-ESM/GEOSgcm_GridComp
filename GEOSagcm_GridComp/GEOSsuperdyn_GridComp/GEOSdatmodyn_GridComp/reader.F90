SUBROUTINE ARM2(FILENAME, NT, NLEVEL, NLAYR,                 &
                PREF_MODEL_E,                                &
                TIME,                                        &
                YY,                                          &
                MO,                                          &
                DD,                                          &
                HH,                                          &
                MM,                                          &
                PCP_OBS,                                     &
                TS_AIR,                                      &
                TG_SOIL,                                     &
                TSKIN,                                       &
                QSFCAIR,                                     &
                QSKIN,                                       &
                PSFC,                                        &
                LHF,                                         & 
                SHF,                                         & 
                ptop,                                        & 
                tt,                                          &
                qq,                                          &
                uu,                                          &
                vv,                                          &
                omega,                                       &
                T_H_adv,                                     &
                T_V_adv,                                     &
                Q_H_adv,                                     &
                Q_V_adv,                                     &           
                useana,                                      &
                T_ana,                                       &
                Q_ana,                                       &
                Q1,                                          &
                Q2,                                          &
                P_MODEL_E                                    )


        use MAPL_Mod
        use GEOS_UtilsMod, only : GEOS_Qsat, GEOS_DQsat

        IMPLICIT NONE
! INPUT        XXXXX
! NT:  time slice number of ARM data
! NLEVEL: vertical levels of ARM data
! NLAYR:  vertical layer of Single column model
! FILENAME: file name of ARM data
!	
        INTEGER,  INTENT(IN) :: Nt, NLEVEL, NLAYR
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        REAL,  DIMENSION(0:NLAYR),     INTENT(IN   ) :: PREF_MODEL_E

! OUTPUT
! single layer
        real, dimension(nt), intent(out) :: time, yy, mo, dd, hh, mm
        real, dimension(nt), intent(out) :: pcp_obs, tskin, qsfcair,&
        qskin, psfc, shf, lhf, ts_air, tg_soil
        integer, intent(in) :: useana
        real, intent(out) :: ptop
! multiple-layer
        real, dimension(nt,nlayr),intent(out):: tt, qq, uu, vv, T_H_adv, T_V_adv
        real, dimension(nt,nlayr),intent(out):: Q_H_adv, Q_V_adv, Q1, Q2

        real, dimension(nt,nlayr),intent(out):: T_ana, Q_ana

        real, dimension(nt,0:nlayr), intent(out):: omega, P_MODEL_E
        real, dimension(nlayr) :: p_model
! temporary arrays
        integer :: IVAR , IVAR_sfc, I, J, K 
        real, allocatable, dimension(:,:,:) :: TMP, DV
        real, allocatable, dimension(:) :: TMP1VAR
        real, allocatable, dimension(:) :: TMP2VAR
        real, allocatable, dimension(:,:) :: TMP_sfc
        real, allocatable, dimension(:,:) :: OMTEMP
        real, allocatable, dimension(:)     :: P
! working variables
        real :: pres, prese, temd, PUPP, PDWN, PRGAT, PUPPK, PDWNK 
        real :: PRESK, TEM, TEMU
        integer :: IGD, IGU, IGTLEV, ITOP, UNIT
        integer :: NT_DATFILE, NLEVEL_DATFILE
        real :: converter

        allocate (P(NLEVEL))
!
        UNIT = GETFILE(FILENAME, form = "formatted")
!
        READ(UNIT,*)
        READ(UNIT,*) NT_DATFILE
        READ(UNIT,*)
        READ(UNIT,*) NLEVEL_DATFILE

        if ((NLEVEL .ne. NLEVEL_DATFILE) .or. (NT .ne. NT_DATFILE)) then
            write(*,*) "***************************************************"
            write(*,*) "SCM NLEVEL or NT in data file doesn't match AGCM.rc "
            write(*,*) "AGCM.rc   NT = ",NT, ", NLEVEL = ",NLEVEL
            write(*,*) "data file NT = ",NT_DATFILE, ", NLEVEL = ",NLEVEL_DATFILE
            write(*,*) "striking in protest"
            stop
        endif

        READ(UNIT,*)
        READ(UNIT,*) ! pressure units
        READ(UNIT,13) P
        READ(UNIT,*)
        READ(UNIT,13) TIME
        READ(UNIT,*)
        READ(UNIT,13) YY
        READ(UNIT,*)
        READ(UNIT,13) MO
        READ(UNIT,*)
        READ(UNIT,13) DD
        READ(UNIT,*)
        READ(UNIT,13) HH
        READ(UNIT,*)
        READ(UNIT,13) MM
        READ(UNIT,*)
        READ(UNIT,*)  IVAR
        
        allocate (TMP(IVAR,NLEVEL,NT))
        allocate (TMP1VAR(NLEVEL))
        allocate (DV (IVAR,NLAYR, NT))
        allocate (TMP2VAR(NLAYR))

        do i = 1, IVAR
          READ(UNIT,*)  ! label
          READ(UNIT,*)  ! unit
          do k = 1, NT
             READ(UNIT,*)  ! time
!             read(UNIT,13) (TMP(I,J,K),K=1,ivar)
             read(UNIT,13) (TMP(I,J,K),j=1,nlevel)
          END DO
        END DO
        READ(UNIT,*)
        READ(UNIT,*) IVAR_sfc
!        READ(UNIT,*) ! temp labels
        allocate (TMP_sfc(IVAR_sfc,NT))
        DO K = 1, ivar_sfc 
           READ(UNIT,*) ! label
           READ(UNIT,*) ! unit
           READ(UNIT,13) (TMP_sfc(K,i),i=1,nt)
        END DO
            write(*,*) "after 2d"
 12     format(5e15.7)
 13     format(e15.7)
 	CALL FREE_FILE(UNIT)

!  Make ihe fields for SCM
!  1.  SURFACE 

        do i = 1, NT
        PCP_OBS(i) = TMP_sfc(5,i)*24.0          ! mm/day
        PSFC(i)    = TMP_sfc(4,i)               ! mb
        TSKIN(i)   = TMP_sfc(3,i)               ! K
        end do

        TG_SOIL    = TSKIN
        TS_AIR     = TSKIN
        do i = 1, NT
        SHF(i)    = TMP_sfc(1,i)
        LHF(I)    = TMP_sfc(2,i)
        end do

        TMP(6,:,:) = TMP(6,:,:)/3600.0          ! HTA k/s
        TMP(7,:,:) = TMP(7,:,:)/3600.0
        TMP(8,:,:) = TMP(8,:,:)/3600.0
        TMP(9,:,:) = TMP(9,:,:)/3600.0

        if ( useana.eq.1  ) then 
        TMP(10,:,:) = TMP(10,:,:)/3600.0 
        TMP(11,:,:) = TMP(11,:,:)/3600.0
        endif
        if (filename.eq."arm_97jul.dat".or.filename.eq."armtwp_ice.dat") then
         TMP(12,:,:) = TMP(12,:,:)/3600.0     ! Horiz s adv (identical to horiz temp adv)
         TMP(13,:,:) = TMP(13,:,:)/3600.0     ! vertical s adv (includes omega*alpha)
        endif
        do i = 1, NT
        do j = 1, NLEVEL
        TMP(1,J,I) = TMP(1,J,I)*((1000.0/P(J))**0.286)
        END DO
        END DO

        do i = 1, NT
        p_model_e(i,0:nlayr) = pref_model_e(0:nlayr)*psfc(i)*100./pref_model_e(nlayr)
        p_model(1:nlayr) = (p_model_e(i,0:nlayr-1) + p_model_e(i,1:nlayr))*0.5/100.
        enddo

!       Interpolate vertically

        do i = 1, NT
        do k = 1, nlayr 
        PRES = P_MODEL(K)
        IF (PRES > P(1) ) THEN
        IGD  = 1
        IGU  = 1
        PUPP = p(1)
        PDWN = PSFC(I)
        GOTO 205
        END IF
        IF (PRES < P(NLEVEL)) THEN
        IGD = NLEVEL
        IGU = NLEVEL
        PUPP = PRES
        PDWN = P(NLEVEL)
        GOTO 205
        END IF
        PUPP = P(NLEVEL)
        PDWN = P(1)
        IGU  = NLEVEL
        IGD  = 2
        DO 204 IGTLEV = 2, NLEVEL
        PRGAT = P(IGTLEV)
        IF (PRGAT <= PRES ) THEN
        PUPP = PRGAT
        IGU  = IGTLEV
        GO TO 900
        END IF
        IF (PRGAT > PRES ) THEN
        PDWN = PRGAT
        IGD  = IGTLEV
        END IF
 204	CONTINUE
 900	CONTINUE
 205	CONTINUE
 	PUPPK = PUPP**0.286
        PDWNK = PDWN**0.286
        PRESK = PRES**0.286

        TEM   = PUPPK - PDWNK
        TEMU  = (PUPPK - PRESK)/TEM
        TEMD  = 1.0 - TEMU
        DO J = 1, IVAR
        DV(J,K,I) = TMP(J,IGU,I)*TEMD + TMP(J,IGD,I)*TEMU
        END DO
        END DO
        END DO

! Need to handle omega differently - interp to model edge pressures

        allocate (OMTEMP (0:NLAYR, NT))

        do i = 1, NT
        do k = 0, nlayr 
        PRES = P_MODEL_E(i,K)/100.
        IF (PRES > P(1) ) THEN
        IGD  = 1
        IGU  = 1
        PUPP = p(1)
        PDWN = PSFC(I)
        GOTO 305
        END IF
        IF (PRES < P(NLEVEL)) THEN
        IGD = NLEVEL
        IGU = NLEVEL
        PUPP = PRES
        PDWN = P(NLEVEL)
        GOTO 305
        END IF
        PUPP = P(NLEVEL)
        PDWN = P(1)
        IGU  = NLEVEL
        IGD  = 2
        DO 304 IGTLEV = 2, NLEVEL
        PRGAT = P(IGTLEV)
        IF (PRGAT <= PRES ) THEN
        PUPP = PRGAT
        IGU  = IGTLEV
        GO TO 800
        END IF
        IF (PRGAT > PRES ) THEN
        PDWN = PRGAT
        IGD  = IGTLEV
        END IF
 304	CONTINUE
 800	CONTINUE
 305	CONTINUE
 	PUPPK = PUPP**0.286
        PDWNK = PDWN**0.286
        PRESK = PRES**0.286

        TEM   = PUPPK - PDWNK
        TEMU  = (PUPPK - PRESK)/TEM
        TEMD  = 1.0 - TEMU
        OMTEMP(k,i) = TMP(5,IGU,I)*TEMD + TMP(5,IGD,I)*TEMU
        END DO
        END DO

         converter = 1000.
        do K = 1, NLAYR
          do I = 1, NT
           tt(i,k) = dv(1,k,I)*((p_model(k)/1000.0)**0.286)
           qq(i,k) = dv(2,K,I)*converter ! kg/kg => g/kg
           uu(i,k) = dv(3,K,I)
           vv(i,k) = dv(4,K,I)
           T_H_ADV(i,k) = -dv(6,K,I)
!=============================================================================
       if (filename.eq."arm_97jul.dat".or.filename.eq."armtwp_ice.dat") then
           T_V_ADV(i,k) = -dv(13,K,I)
        elseif (filename.eq."arm_kwjx.dat"  .or. filename.eq."TRMM_LBA.dat") then
           T_V_ADV(i,k) = -dv(7,K,I) - MAPL_RGAS * tt(i,k) * omtemp(k,i)/3600. / (MAPL_CP * p_model(k) )
        else
           T_V_ADV(i,k) = -dv(7,K,I)
        endif

           Q_H_ADV(i,k) = -dv(8,K,I)
           Q_V_ADV(i,k) = -dv(9,K,I)

        if ( useana.eq.1  ) then
           T_ana(i,k) = -dv(10,K,I)   
           Q_ana(i,k) = -dv(11,K,I)
        endif
        END DO
        END DO

        do K = 0, NLAYR
          do I = 1, NT
           omega(i,k) = omtemp(K,I) / 36.  ! convert from mb/hour to Pa/sec
           END DO
        END DO

        Q1 = -999.0
        Q2 = -999.0

        do i = 1, NT
        ITOP = 1
        do k = 1, nlayr-1
        if((p_model(k)<=minval(p)) .AND. (p_model(k+1)>minval(p))) ITOP = K
        END DO
        do k = 1, ITOP
        tt(i,k) = tt(i,itop)
        end do
        end do

        do i = 1, NT
        QSKIN(I) = 1000.*GEOS_QSAT (TSKIN(i),PSFC(i))
        QSFCAIR(i) = QQ(i,nlayr)
        END DO

        ptop = minval(p)*100.      !  top data pressure in Pa

        END SUBROUTINE ARM2
