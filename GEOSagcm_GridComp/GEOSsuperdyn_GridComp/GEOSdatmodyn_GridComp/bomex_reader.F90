SUBROUTINE BOMEX(NT, NLEVEL, NLAYR,                           &
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
                 tt,                                          &
                 qq,                                          &
                 uu,                                          &
                 vv,                                          &
                 T_H_adv,                                     &
                 T_V_adv,                                     &
                 Q_H_adv,                                     &
                 Q_V_adv,                                     &           
                 Q1,                                          &
                 Q2,                                          &
                 P_MODEL_E                                    )

        use MAPL_Mod
        use GEOS_UtilsMod, only : GEOS_Qsat, GEOS_DQsat

        IMPLICIT NONE
! INPUT        XXXXX
! NT:  time slice number of data
! NLEVEL: vertical levels of data
! NLAYR:  vertical layer of Single column model
!	
	INTEGER,  INTENT(IN) :: Nt, NLEVEL, NLAYR
	REAL, DIMENSION(0:NLAYR), INTENT(IN) :: PREF_MODEL_E

! OUTPUT
! single layer
	real, dimension(nt), intent(out) :: time, yy, mo, dd, hh, mm
	real, dimension(nt), intent(out) :: pcp_obs, tskin, qsfcair,&
	qskin, psfc, shf, lhf, ts_air, tg_soil
! multiple-layer
	real, dimension(nt,nlayr),intent(out):: tt, qq, uu, vv, &
	T_H_adv, T_V_adv,    &
        Q_H_adv, Q_V_adv, Q1, Q2
	real, dimension(nt,0:nlayr), intent(out):: P_MODEL_E
	real, dimension(nlayr) :: p_model
! temporiy array
	integer :: IVAR_sfc, I, J, K 
	integer, parameter :: IVAR = 9
	real :: SHF_C, LHF_C
	real, dimension(9,NLEVEL,NT) :: TMP,DV
! working variables
	real :: pres, temd, PUPP, PDWN, PRGAT, PUPPK, PDWNK 
	real :: pmass, PRESK, TEM, TEMU, WK1_C
	integer :: IGD, IGU, IGTLEV, ITOP, UNIT
	real, dimension(75) :: P, Z
! 
        data Z/                                                         &
             20.00,      60.00,     100.00,     140.00,     180.00,     &
            220.00,     260.00,     300.00,     340.00,     380.00,     &
            420.00,     460.00,     500.00,     540.00,     580.00,     &
            620.00,     660.00,     700.00,     740.00,     780.00,     &
            820.00,     860.00,     900.00,     940.00,     980.00,     &
           1020.00,    1060.00,    1100.00,    1140.00,    1180.00,     &
           1220.00,    1260.00,    1300.00,    1340.00,    1380.00,     &
           1420.00,    1460.00,    1500.00,    1540.00,    1580.00,     &
           1620.00,    1660.00,    1700.00,    1740.00,    1780.00,     &
           1820.00,    1860.00,    1900.00,    1940.00,    1980.00,     &
           2020.00,    2060.00,    2100.00,    2140.00,    2180.00,     &
           2220.00,    2260.00,    2300.00,    2340.00,    2380.00,     &
           2420.00,    2460.00,    2500.00,    2540.00,    2580.00,     &
           2620.00,    2660.00,    2700.00,    2740.00,    2780.00,     &
           2820.00,    2860.00,    2900.00,    2940.00,    2980.00/        
        data P/                                                         &
         101271.27,  100814.89,  100359.99,   99906.54,   99454.52,     &
          99003.97,   98554.85,   98107.19,   97660.96,   97216.16,     &
          96772.80,   96330.87,   95890.38,   95451.36,   95013.89,     &
          94578.03,   94143.77,   93711.06,   93279.96,   92850.43,     &
          92422.48,   91996.09,   91571.27,   91148.02,   90726.32,     &
          90306.16,   89887.56,   89470.50,   89054.99,   88641.00,     &
          88228.54,   87817.62,   87408.20,   87000.30,   86593.91,     &
          86189.03,   85785.66,   85383.85,   84983.80,   84585.54,     &
          84189.09,   83794.42,   83401.54,   83010.41,   82621.05,     &
          82233.46,   81847.59,   81463.48,   81081.08,   80700.41,     &
          80321.39,   79943.88,   79567.80,   79193.13,   78819.89,     &
          78448.09,   78077.70,   77708.70,   77341.12,   76974.96,     &
          76610.20,   76246.81,   75884.82,   75524.23,   75165.02,     &
          74807.18,   74450.72,   74095.63,   73741.90,   73389.53,     &
          73038.53,   72688.87,   72340.56,   71993.60,   71647.98/      
	  P(:) = 0.01*P(:)
	TIME = -999.0
	YY   = 1969.0
	MO   = 6.0
	MM   = 0.0
	do i = 1, 5
	do j = 1, 24
	k = j + (i-1)*24
	dd(k) = 21 + i
	hh(k) = j - 1
	end do
	end do
!  1.  SURFACE 
	
        do i = 1, NT
	PCP_OBS(i) = -999.0                     ! mm/day
	PSFC(i)    = 1015.00                    ! mb
	TSKIN(i)   = 300.375                    ! K
	end do
	TG_SOIL    = TSKIN
	TS_AIR     = TSKIN
	SHF_C = 1004.0*PSFC(1)*100.0/(287.0*TSKIN(1))
	LHF_C =  2.5E6*PSFC(1)*100.0/(287.0*TSKIN(1))
	do i = 1, NT
	SHF(i)    = SHF_C*8.0E-3
	LHF(I)    = LHF_C*5.2E-5
	end do
! Tendency due to advection
	TMP(6,:,:) = 0.0         ! H-Temp. Advection 
	TMP(7,:,:) = 0.0 
	TMP(9,:,:) = 0.0         ! V-Q. Advection

	do i = 1, NT
	do j = 1, NLEVEL
	if(Z(j)>=0.0.and.Z(j)<300.0) TMP(8,j,I) = -1.2*1.0E-8
	if(Z(j)>=300.0.and.Z(j)<500.0) &
	TMP(8,j,I) = -(1.2E-8-1.2E-8*(Z(j)-300.0)/200.0)
	if(Z(j) >= 500.0) TMP(8,j,I) = 0.0
	end do
	end do
	TMP(8,:,:) = -TMP(8,:,:)*1.E+3   ! converting into g/kg/s
! Temperature
	do i = 1, NT
	do j = 1, NLEVEL
	if(Z(j)>=0.0.and.Z(j)<520.0) TMP(1,J,I) = 298.7
	if(Z(j)>=520.0.and.Z(j)<1480.0)  &
	TMP(1,J,I) = 298.7+(302.4-298.7)/(1480-520)*(z(j)-520.0)
	if(Z(j)>=1480.0.and.Z(j)<2000.0) &
	TMP(1,J,I) =302.4+(308.2-302.4)/(2000-1480)*(z(j)-1480.0)
	if(Z(j)>=2000.0) TMP(1,J,I)=308.2+3.65E-3*(z(j)-2000.0)
	END DO
	END DO
	do i = 1, NT
	do j = 1, NLEVEL
	TMP(1,J,I) = TMP(1,J,I)*((P(J)/1000.0)**0.286)
	end do
	end do
! Moisture
	do i = 1, NT
	do j = 1, NLEVEL
	if(Z(j)>=0.0.and.Z(j)<520.0)&
	TMP(2,J,I) = 17.0 + (16.3 - 17.0)/(520)*Z(J)
	if(Z(j)>=520.0.and.Z(j)<1480.0)  &
	TMP(2,J,I) = 16.3+(10.7-16.3)/(1480-520)*(z(j)-520) 
	if(Z(j)>=1480.0.and.Z(j)<2000.0) &
	TMP(2,J,I) = 10.7+(4.2-10.7)/(2000-1480)*(z(j)-1480)
	if(Z(j)>=2000.0) TMP(2,J,I) = 4.2-1.2E-3*(z(j)-2000)
	end do
	end do
! U-Wind
	do i = 1, NT
	do J = 1, NLEVEL
	if(Z(j)>=0.0.and.Z(j)<700.0)&
	TMP(3,J,I) = -8.75
	if(Z(j)>=700.0) &
	TMP(3,J,I) = -8.75 + 1.8E-3*(z(j)-700)
	end do
	end do
! V-Wind
	TMP(4,:,:) = 0.0
! Omega-Wind
	! W (dZ/dt, m/s)
	do i = 1, NT
	do J = 1, NLEVEL
	if(Z(j)>=0.0.and.Z(j)<1500.0)&
	TMP(5,J,I) = - (0.0065/1500)*z(j)
	if(Z(j)>=1500.0.and.Z(j)<2100.0) &
	TMP(5,J,I) = -0.0065+0.0065/(2100-1500)*(z(j)-1500)
	if(Z(j)>=2100.0 ) TMP(5,J,I) = 0.0
	end do
	end do
! Adiabtic term in Temp. Vectical Advection term
	do I = 1, NT
	do J = 1, NLEVEL
	TMP(7,J,I)=-9.8*TMP(5,J,I)/1004.0
	end do
	end do

	! convertint into Omega
	do i = 1, NT
	do J = 1, NLEVEL
	WK1_C = -P(j)*980.0/(287.0*TMP(2,J,I))
	TMP(5,J,I) = WK1_C*TMP(5,J,I)
	end do
	end do

!       Interpolate vertically
	do i = 1, NT
	p_model_e(i,0:nlayr) = pref_model_e(0:nlayr)*100.*psfc(i)/pref_model_e(nlayr)
	p_model(1:nlayr) = (p_model_e(i,0:nlayr-1) + p_model_e(i,1:nlayr))*0.5/100.
	
	do k = 1, nlayr 
	PRES = P_MODEL(K)
	IF (PRES > P(1) ) THEN
	IGD  = 1
	IGU  = 1
	PUPP = 965.0
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
	do K = 1, NLAYR
	   do I = 1, NT
	   tt(i,k) = dv(1,k,I)*((1000.0/p_model(k))**0.286)
	   qq(i,k) = dv(2,K,I)
	   uu(i,k) = dv(3,K,I)
	   vv(i,k) = dv(4,K,I)
	   T_H_ADV(i,k) = -dv(6,K,I)
	   T_V_ADV(i,k) = -dv(7,K,I)
	   Q_H_ADV(i,k) = -dv(8,K,I)
	   Q_V_ADV(i,k) = -dv(9,K,I)
	   END DO
        END DO
	Q1 = -999.0
	Q2 = -999.0
	do i = 1, NT
	ITOP = 1
	do k = 1, nlayr
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
	END SUBROUTINE BOMEX
