      SUBROUTINE STRATUS_IC(                                &
                     nlayr,                              &
                     PREF_MODEL_E,                       &
                     STRATUS_SST,                        &
                     STRATUS_OMG,                        &
                     STRATUS_SFCAIR_TEMP,                &
                     time,                               &
                     yy,                                 &
                     mo,                                 &
                     dd,                                 &
                     hh,                                 &
                     mm,                                 & 
                     PCP_OBS,                            &
                     TS_AIR,                             &
                     TG_SOIL,                            &
                     TSKIN,                              &
                     QSFCAIR,                            &
                     QSKIN,                              &
                     PSFC,                               &
                     tt,                                 &
                     qq,                                 &
                     uu,                                 &
                     vv,                                 &
                     T_H_adv,                            &
                     T_V_adv,                            &
                     Q_H_adv,                            &
                     Q_V_adv,                            &
                     Q1,                                 &
                     Q2,                                 & 
                     OMG,                                & 
                     P_MODEL_E                           )

        use GEOS_Mod

        IMPLICIT None


! ARGUMENTS
! INPUT

      INTEGER, parameter                           :: NT=1

      INTEGER,                       INTENT(IN   ) :: nlayr
      REAL,  DIMENSION(0:NLAYR),     INTENT(IN   ) :: PREF_MODEL_E
      REAL,  INTENT(IN   ) :: STRATUS_SFCAIR_TEMP, STRATUS_SST , STRATUS_OMG 

! OUTPUT
! Singli-layer:
      real, DIMENSION(nt       ), INTENT(  OUT) :: time       !Calenday day
      real, DIMENSION(nt       ), INTENT(  OUT) ::  yy        !Year
      real, DIMENSION(nt       ), INTENT(  OUT) ::  mo        !Month
      real, DIMENSION(nt       ), INTENT(  OUT) ::  dd        !Day
      real, DIMENSION(nt       ), INTENT(  OUT) ::  hh        !Hour
      real, DIMENSION(nt       ), INTENT(  OUT) ::  mm        !Minutes
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: PCP_OBS
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: TS_AIR 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: TG_SOIL 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: TSKIN 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: QSKIN 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: QSFCAIR 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: PSFC

! Multiple-layer

      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: tt
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: qq
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: uu
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: vv
!     REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: ww  
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: T_H_adv
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: T_V_adv
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: Q_H_adv
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: Q_V_adv
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: Q1
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: Q2
      REAL, DIMENSION(NT, 0:nlayr), INTENT(INOUT) :: OMG
      REAL, DIMENSION(NT, 0:nlayr), INTENT(  OUT) :: P_MODEL_E

      REAL, DIMENSION(NLAYR)       :: THETA00, TEMP, Q, PLO, PKO, DZ
      REAL, DIMENSION(0:NLAYR)     :: ZE

      real :: r1,r2,r3,r4,r5,r6 ,lat0,lon0,hour,omeg0
      integer :: i0,i1,L,LM,iyyyy,imo,idd,ihh,imm,iss,I,UNIT,n


       P_MODEL_E(1,:) = PREF_MODEL_E(:)
       PLO   = ( P_MODEL_E(1,1:NLAYR) +  P_MODEL_E(1,0:NLAYR-1) )/2.0
       PKO   = ( PLO / MAPL_P00 )**MAPL_KAPPA
 
       THETA00= STRATUS_SFCAIR_TEMP

       THETA00 = THETA00 + (100000.- PLO)*30./90000.

       TEMP   = THETA00*PKO

       where( TEMP<190.)
         TEMP=190.
       endwhere

       Q      = 0.8*GEOS_QSAT (TEMP, PLO/100.)


#ifdef POLAR 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    ! attempt to get polar khrad prob
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       where (PLO > 85000. )
         THETA00=250.
       endwhere

       where ( PLO <= 85000. )
         THETA00= 260.*( (85000./plo)**(1./9.) )
       endwhere

       where( THETA00 > 1200. )
         THETA00=1200.
       endwhere

       TEMP   = THETA00*PKO

       Q      = 0.8*GEOS_QSAT (TEMP, PLO/100.)
       where( (PLO <=185000.) .and. (PLO > 85000.) )
          Q      = 1.15*GEOS_QSAT (TEMP, PLO/100.)
       endwhere
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
#endif

       TT(1,:)= TEMP
       QQ(1,:)= Q
       

       iyyyy =  1900
       imo   =  7
       idd   =  1



       yy    =  1.0*iyyyy
       mo    =  1.0*imo   
       dd    =  1.0*idd
       hh    =  0.00
       mm    =  0.00
      

       time       =  -999.9
       PCP_OBS    = -999.9
       TS_AIR     = THETA00*( ( P_MODEL_E(1,NLAYR) / MAPL_P00 )**MAPL_KAPPA )
       TG_SOIL    = -999.9 
       PSFC       = P_MODEL_E(1,NLAYR)
       TSKIN      = STRATUS_SST
       QSFCAIR    = QQ(1,NLAYR)
       QSKIN      = GEOS_QSAT (TSKIN , PSFC/100.)

      uu = 0.0 
      vv = 0.0
      T_H_adv =0.
      T_V_adv =0.
      Q_H_adv =0.
      Q_V_adv =0.
      Q1      =0.
      Q2      = 0.
      OMEG0   = STRATUS_OMG * 100./86400.

      do n=1,nt
           omg(n,:) = ( p_model_e(n,:)-psfc(n))*( 10000.-p_model_e(n,:) ) &
                       /( ( (psfc(n)-10000.)/2)**2 )
      end do

      do n=1,nt
           omg(n,:) = ( p_model_e(n,:)-psfc(n))*( 60000.-p_model_e(n,:) ) &
                       /( ( (psfc(n)-60000.)/2)**2 )
      end do

      where(omg < 0. )
        omg=0.
      endwhere

      omg = OMEG0 *  omg 

      do n=1,nlayr
         vv(:,n) = ((plo(n) - psfc(:) )/10000.)**2
      enddo

      vv=min(vv,5.0) +3.0
      uu=5.0


      uu = 1.
      vv = 1.

	write(*,*) 'Done!'
      END SUBROUTINE STRATUS_IC
