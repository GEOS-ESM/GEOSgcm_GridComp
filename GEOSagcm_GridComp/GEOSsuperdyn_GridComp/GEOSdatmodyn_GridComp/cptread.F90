SUBROUTINE CPTREAD(FILENAME, NT, NLEVEL, NLAYR,                 &
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
                OMEGA,                                       &
                P_MODEL_E                                    )

        use GEOS_Mod

        IMPLICIT NONE
! INPUT        XXXXX
! NT:  time slice number of data
! NLEVEL: vertical levels of data
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
! multiple-layer
	real, dimension(nt,nlayr),intent(out):: tt, qq, uu, vv, &
	T_H_adv, T_V_adv,    &
        Q_H_adv, Q_V_adv
	real, dimension(nt,0:nlayr), intent(out):: P_MODEL_E, OMEGA
	real, dimension(nlayr) :: p_model
! temporiy array
	integer :: IVAR , IVAR_sfc, I, J, K , L, NVAR, L2
	real, dimension(nlevel)          :: P_data
	real, dimension(nlayr,nt)        :: DV
	real, dimension(nt,nlayr)        :: OMEGL
	real, dimension(nlevel,nt)       :: TMP
	real, parameter :: ptop = 10.

	real :: pres, temd, PUPP, PDWN, PRGAT, PUPPK, PDWNK 
	real :: pmass, PRESK, TEM, TEMU
	integer :: IGD, IGU, IGTLEV, ITOP, UNIT, NLEVSFILE, NTFILE

!
	UNIT = GETFILE(FILENAME, form = "unformatted")
!

	READ(UNIT) YY
          write(*,*) 'yy'
	READ(UNIT) MO
          write(*,*) 'mo'
	READ(UNIT) DD
          write(*,*) 'dd'
	READ(UNIT) HH
          write(*,*) 'hh'
	READ(UNIT) MM
          write(*,*) 'mm'
	READ(UNIT) P_data
          write(*,*) 'pdata'
	READ(UNIT) PSFC
          write(*,*) 'psfc'
	READ(UNIT) TSKIN
          write(*,*) 'tskin'

	do i = 1, NT
           p_model_e(i,0:nlayr) = pref_model_e(0:nlayr) * psfc(i) / pref_model_e(nlayr)
        end do
   	
        p_model  =  (p_model_e(1,0:nlayr-1) + p_model_e(1,1:nlayr))*0.5 
                  write(212) p_model
        IVAR = 7

	do nvar = 1, IVAR
	   read(UNIT) tmp
              write(*,*) 'ivar=',nvar

           do i = 1, NT
  	      p_model  =  (p_model_e(i,0:nlayr-1) + p_model_e(i,1:nlayr))*0.5 
              !! slop     = ( TMP(nlevel-1,i) - TMP(nlevel,i) ) / ( p_data(nlevel-1) - p_data(nlevel) )
            
              do l2 = 1,nlayr

                 if( p_model(l2) > p_data(nlevel) ) then
                 ! where( p_model > p_data(nlevel) )
                   dv( l2 , i ) = tmp(nlevel,i)   !+ slop*( p_data(nlevel) - p_model )
                 ! endwhere
                 endif

                 do l = 1,nlevel-1
                    if(  ( p_model(l2) < p_data(l+1) ) .AND. ( p_model(l2) >= p_data(l) ) ) then
                    !!where(  ( p_model < p_data(l+1) ) .AND. ( p_model >= p_data(l) ) )
                         dv( l2 , i ) = tmp(l,i) + (p_model(l2)-p_data(l))*(tmp(l+1,i)-tmp(l,i)) / (p_data(l+1)-p_data(l))
                    !!endwhere
                    endif
                 end do

                 if( p_model(l2) <= p_data(1) ) then
                 !!where( p_model <= p_data(1) )
                     dv( l2 , i ) = tmp(1,i)  
                 !!endif
                 endif

              end do  ! over L2

            end do ! over i


                   write( 212 ) DV

           SELECT CASE( NVAR )

                CASE( 1 )
                    TT = TRANSPOSE( DV )                    
                CASE( 2 )
                    QQ = TRANSPOSE( DV ) * 1000.                    
                CASE( 3 )
                    UU = TRANSPOSE( DV )                    
                CASE( 4 )
                    VV = TRANSPOSE( DV )                    
                CASE( 5 )
                    OMEGL   = TRANSPOSE( DV )                    
                CASE( 6 )
                    T_H_ADV = TRANSPOSE( DV )                    
                CASE( 7 )
                    Q_H_ADV = TRANSPOSE( DV ) * 1000.                   

           END SELECT

	END DO


	TG_SOIL    = TSKIN
	TS_AIR     = TSKIN
	QSKIN      = 1000.*GEOS_QSAT (TSKIN ,PSFC/100. )
	QSFCAIR    = QQ(:,nlayr)
        PCP_OBS    = -99999.
        SHF        = -99999.
        LHF        = -99999.   
        TIME       = -99999.	


        DO L=1,nlevel-1
           OMEGA(:,L) = 0.5*( OMEGL(:,L) + OMEGL(:,L+1) )
        END DO
        OMEGA(:,0)=0.0
        OMEGA(:,NLEVEL)=0.0


	END SUBROUTINE CPTREAD
