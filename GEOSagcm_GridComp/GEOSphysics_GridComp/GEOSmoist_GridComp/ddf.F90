! $Id$
MODULE DDF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This will be a simple downdraft routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use ESMF
use MAPL_Mod
use GEOS_UtilsMod, only : DQSAT=>GEOS_DQsat, QSAT=>GEOS_Qsat


 IMPLICIT NONE

 PRIVATE

 PUBLIC T_prec_partition
 PUBLIC T_DDF_CTL
 PUBLIC DDF1

type T_prec_partition
   character*72 :: type
   real         :: ratio
   real         :: min_ratio
   real         :: max_ratio
end type T_prec_partition

type T_DDF_CTL
   real    :: ALPHA_DDF_UDF
   real    :: AREAL_FRACTION
   logical :: THETA_IS_THVAR
   type(T_prec_partition) :: partition
end type T_DDF_CTL

real, parameter :: RHO_W  =  1.0e3  ! Density of liquid water in kg/m^3

 CONTAINS
!========================================================================

  SUBROUTINE DDF1( CTLS , DT, T , QV, QL, QI, QPL, QPI , Z, ZE, P, PE, PK , UMFC, & 
                   ZSCALE_DIAG, DQDT_DIAG, DTDT_DIAG , DMFC_DIAG )

    type( T_DDF_CTL ), intent(inout)       :: CTLS   

    real  ,  intent(in)                    :: DT
    real  ,  intent(in)    , dimension(0:) :: PE, ZE, UMFC
    real  ,  intent(in)    , dimension(:)  :: Z, P, PK
    real  ,  intent(inout) , dimension(:)  :: T, QV, QL, QI, QPL, QPI

    real  ,  pointer,        dimension(:)  :: DQDT_DIAG, DTDT_DIAG, DMFC_DIAG
    real                                   :: ZSCALE_DIAG
    
    real  , dimension(0:size(T,1) ) :: DMFC,ZEr,WE

    real  , dimension(1:size(T,1) ) :: DM,DZ, SS, DDFRAC, rho
    real  , dimension(1:size(T,1) ) :: Tc, QVc, QLc, QIc, QPLc, QPIc, SSc, QSc
    real  , dimension(1:size(T,1) ) :: RHc, DQSc
    real  , dimension(1:size(T,1) ) :: EVAP_DIAG, SUBL_DIAG, RH_DIAG

    real  :: Zscale , DMFCx, QVcL, AA, BB, Efactor,FlakRad,DropRad,FallSpLiq,FallSpIce,TinLayer,Evap,Subl,Efactor2

    real  :: TVE0, TVE1, TVQ0, TVQ1 , CFLchk,zdet

    integer :: LM, L, Ltop
  
    logical, parameter :: POST_CORRECTION =.TRUE.
    logical, parameter :: BUDGET_CHECKS   =.TRUE.

! Begin 
!---------

    LM = size(T,1)
                                       ! layer mass (kg m-2)
    DM  = ( PE(1:LM) -  PE(0:LM-1) )*100./MAPL_GRAV 

    DZ  = -( ZE(1:LM) -  ZE(0:LM-1) )  ! layer thickness (m)

    rho = DM/DZ                        ! layer density (kg m-3)
      
    ZEr = ZE - ZE(LM)                  ! relative edge geometric heights (m)
 
    if(CTLS%THETA_IS_THVAR)  T = PK * T 

      ! normalized downdraft profile DMF / MAX(DMF)
      !--------------------------------------------
   Zscale = fzscale( CTLS, ZEr, UMFC )   ! 5000.

    where ( ZEr < Zscale )
         DMFC = ( ZEr / Zscale ) * ( (1.0 - ZEr/Zscale)**2 ) 
    elsewhere
         DMFC = 0.0
    endwhere
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! <<flattened>> downdraft
        !----------------------------------
!    zscale = max( zscale, 3000. )
!
!    zdet=min( zscale*0.25, 1000. )
!
!    DMFC = 0.0
!    where( ZEr < Zdet )
!        DMFC = ZEr / Zdet
!    endwhere
!
!    where( ( ZEr > Zscale-Zdet) .and. (ZEr <= Zscale ) )
!        DMFC = (Zscale - Zer ) / Zdet 
!    endwhere
!
!    where( ( ZEr <= Zscale-Zdet) .and. (ZEr > Zdet ) )
!        DMFC = 1.000
!    endwhere
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if(MAXVAL(DMFC)>0.0) DMFC = DMFC / MAXVAL( DMFC )

      ! find mass flux scaling
      !-------------------------
    call closure_1( UMFC , DMFCx , CTLS  )
    DMFC=-DMFC*DMFCx
   

            ! CFL check  (hardwired at 0.5)
            !--------------------------------------
           CFLchk = MAXVAL( -DMFC(0:LM-1) * DT / DM )
           if (CFLchk > 0.5) DMFC = ( 0.5 / CFLchk )*DMFC



!----------------------------------------------------
! At this point correctly-scaled downdraft mass flux
! profile should exist with units of kg m-2 s-1 
!----------------------------------------------------      
                  if (associated(DMFC_DIAG))   DMFC_DIAG   = DMFC
                                               zscale_diag = zscale



      ! calculate and save mass wt vertical integrals
      ! for later post facto corrections
      !-------------------------------------------
    TVQ0    = SUM( (  QV +  QPL + QPI  )*DM , 1 )
    TVE0    = SUM( (  MAPL_CP*T + MAPL_ALHL*QV - MAPL_ALHF*QPI ) *DM , 1 )


      ! remove some precip condensate for draft
      !------------------------------------------
    call precipfraction( CTLS, Zscale, ZEr, QPL, QPLc )
    call precipfraction( CTLS, Zscale, ZEr, QPI, QPIc )
 
      ! Now <<amplify>> QPs within rain shaft
      !------------------------------------------
    QPLc     = QPLc / CTLS%AREAL_FRACTION
    QPIc     = QPIc / CTLS%AREAL_FRACTION
      

      
      ! make up a vertical velocity profile by assuming
      ! a constant cross-sectional area
      !---------------------------------------------    
    WE(1:LM)  = DMFC(1:LM) / ( rho * CTLS%AREAL_FRACTION )


!-----------------------------------------------------
! Now find in-draft profiles
!-----------------------------------------------------

         !* first locate top of downdraft 
         !--------------------------------
       Ltop=1
       do L=1,LM
          if ( ABS( DMFC(L) ) > 0. ) then
             Ltop = L
             exit
          endif
       enddo       

       Ltop=6

      ! Now start calc-ing from Ltop down
      !-----------------------------------
    QVc = QV
    SS  = MAPL_CP*T + MAPL_GRAV*Z
    SSc = SS
    Tc  = T
    do L=Ltop,LM

         QVc(L) = fx1( QVc(L-1:L) , QV(L), DZ(L) , DMFC(L-1:L) )
         SSc(L) = fx1( SSc(L-1:L) , SS(L), DZ(L) , DMFC(L-1:L) )

         Tc(L)  = ( SSc(L) - MAPL_GRAV*Z(L) )/ MAPL_CP
          
         DQSc(L) = DQSAT( Tc(L) , P(L), qsat=QSc(L) )
         RHc(L)  = QVc(L)/QSc(L)
            

              !-----------------------
                  RH_DIAG(L) = RHc(L)
              !-----------------------

         call MICRO_AA_BB(Tc(L),P(L),QSc(L),AA,BB)

         Efactor = ( 1.00 - RHc(L) ) / ( RHO_W * ( AA  + BB  ))  ! / (1.00 - RHc(L) )
         Efactor = MAX( Efactor, 0.0 )

             ! add precip from above
             !---------------------------------------
         QPLc(L)   = QPLc(L) + DM(L-1)*QPLc(L-1)/DM(L) 
         QPLc(L-1) = 0.00

         QPIc(L)   = QPIc(L) + DM(L-1)*QPIc(L-1)/DM(L) 
         QPIc(L-1) = 0.00

         FallSpLiq = 5.0 ! m/s
         FallSpIce = 5.0 ! m/s

         DropRad   = 0.0002 ! m
         FlakRad   = 0.0002 ! m


         
              ! Evaporate some rain into ddf
         TinLayer =   DZ(L)/FallSpLiq
         Efactor2 =   Efactor / ( 1.0 + DQSc(L)*MAPL_ALHL/MAPL_CP )/ (DropRad**2)
         EVAP     =   QPLc(L)*(1.0 - EXP( -TinLAYER*Efactor2  ) )
              ! Sublimate some hail/snow into ddf
         TinLayer =   DZ(L)/FallSpIce
         Efactor2 =   Efactor / ( 1.0 + DQSc(L)*MAPL_ALHS/MAPL_CP )/ (FlakRad**2)
         SUBL     =   QPIc(L)*(1.0 - EXP( -TinLAYER*Efactor2  ) )

              ! Add effects of phase changes into in-cloud ddf quantities
         QPLc(L)  = QPLc(L) - EVAP
         QPIc(L)  = QPIc(L) - SUBL
         QVc(L)   = QVc(L)  + EVAP + SUBL
         Tc(L)    = Tc(L)   - EVAP * MAPL_ALHL / MAPL_CP - SUBL * MAPL_ALHS / MAPL_CP
              ! Remember to recalculate in-cloud static Energy
         SSc(L)   = MAPL_CP*Tc(L) + MAPL_GRAV*Z(L)

     
              !-----------------------
                  EVAP_DIAG(L) = EVAP
                  SUBL_DIAG(L) = SUBL
              !-----------------------


    end do

               if(BUDGET_CHECKS) TVQ1    = SUM( (  QV +  QPL + QPI  )*DM , 1 )
               if(BUDGET_CHECKS) TVE1    = SUM( (  MAPL_CP*T + MAPL_ALHL*QV  - MAPL_ALHF*QPI ) *DM , 1 )

      ! Return draft precip condensate to grid
      ! box (alternative make sfc precip here?)
      !------------------------------------------
    QPL = QPL + QPLc* CTLS%AREAL_FRACTION
    QPI = QPI + QPIc* CTLS%AREAL_FRACTION

               if(BUDGET_CHECKS) TVQ1    = SUM( (  QV +  QPL + QPI  )*DM , 1 )
               if(BUDGET_CHECKS) TVE1    = SUM( (  MAPL_CP*T + MAPL_ALHL*QV - MAPL_ALHF*QPI ) *DM , 1 )

!--------------------------------------------------------
! At this point we have in-cloud profiles for everything.
! Now calculate env tendencies due to ddf mass flux
!---------------------------------------------------------



              !-----------------------
                  if (associated(DQDT_DIAG)) DQDT_DIAG  = Qv
                  if (associated(DTDT_DIAG)) DTDT_DIAG  = T 
              !-----------------------
    
    do L=Ltop,LM

         QV(L) = eup1( QV(L-1:L) , QVc(L), DZ(L), RHO(L), DT , DMFC(L-1:L) )
         SS(L) = eup1( SS(L-1:L) , SSc(L), DZ(L), RHO(L), DT , DMFC(L-1:L) )

         T(L)  = ( SS(L) - MAPL_GRAV*Z(L) )/ MAPL_CP

    enddo
              !-----------------------
                  if (associated(DQDT_DIAG)) DQDT_DIAG  = ( Qv - DQDT_DIAG )/DT
                  if (associated(DQDT_DIAG)) DTDT_DIAG  = ( T  - DTDT_DIAG )/DT
              !-----------------------


             ! proportional mass corrections
             !  - adds tracer mass proportionally to mixing ratio
             !    to minimize spurious vertical redistributions 
             !  - correct QV first then Static Energy with 
             !    corrected QV
             !_________________________________________________
    if(POST_CORRECTION) then
      TVQ1    = SUM( (  QV +  QPL + QPI )*DM , 1 )
      QV  = QV - (TVQ1 - TVQ0)*QV/TVQ1

      TVE1    = SUM( (  MAPL_CP*T + MAPL_ALHL*QV - MAPL_ALHF*QPI ) *DM , 1 )
      T  = T   - (TVE1 - TVE0)*T / TVE1
    endif

              if(BUDGET_CHECKS) TVE1    = SUM( (  MAPL_CP*T + MAPL_ALHL*QV - MAPL_ALHF*QPI ) *DM , 1 )
              if(BUDGET_CHECKS) TVQ1    = SUM( (  QV +  QPL + QPI  )*DM , 1 )
 

    if(CTLS%THETA_IS_THVAR)  T = T / PK 

  end SUBROUTINE DDF1 


!----------------------------------------------------------------------
!-----------------------------------------------------------------------


  SUBROUTINE precipfraction ( CTLS, Zscale, ZEr, QP, QPc )

    real  ,  intent(in)                      :: Zscale
    real  ,  intent(in)       , dimension(:) :: ZEr(0:)
    real  ,  intent(inout)    , dimension(:) :: QP, QPc
    type(T_DDF_CTL), intent(in)              :: CTLS
 
    real, dimension( size(QP) ) :: ddfrac   

    ddfrac = 0.0
 
    if (TRIM( CTLS%PARTITION%TYPE ) == "UNIFORM" ) then
       ddfrac = CTLS%PARTITION%RATIO 
    endif
    if (TRIM( CTLS%PARTITION%TYPE ) == "HEIGHT_DEP" ) then
       where ( ZEr(1:) > 0.2 *Zscale )
          ddfrac = CTLS%PARTITION%MAX_RATIO
       elsewhere 
          ddfrac = CTLS%PARTITION%MIN_RATIO  + &
                     ZEr(1:)*( CTLS%PARTITION%MAX_RATIO - CTLS%PARTITION%MIN_RATIO ) & 
                           / ( 0.2*Zscale )
       endwhere
    endif

    QPc     = QP * DDFRAC
    QP      = QP - QPc

 
  end SUBROUTINE precipfraction 
!----------------------------------------------------------------------
!-----------------------------------------------------------------------


  SUBROUTINE CLOSURE_1 ( UMFC, DMF0 , CTLS )

    real  ,  intent(in)    , dimension(0:) :: UMFC
    type(T_DDF_CTL), intent(in)            :: CTLS
    real  ,  intent(out)                   :: DMF0


    
    dmf0 = CTLS%ALPHA_DDF_UDF   *MAXVAL( UMFC )

  
  end SUBROUTINE CLOSURE_1 


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

real function fzscale ( ctls, ze, umfc  ) result( zscale ) 
!===================================================
! returns scaling height for downdraft mass flux
! profile (m)
!===================================================

real, intent(in),     dimension(0:)  :: umfc,ze
type(T_DDF_CTL), intent(in)       :: CTLS

real, dimension( 0:size(ze)-1 )   :: ff
real :: deno


ff   =    umfc
!if ( maxval(umfc) > 0. ) then
!   where( umfc > 0.1*maxval(umfc) )
!     ff=1.0
!   endwhere
!endif


deno =    sum( ff , 1 )

if ( deno > 0. ) then
   zscale = SQRT(  sum( ze * ze * ff ,1 ) / deno )
else
   zscale = -1.0
endif

if ( zscale > 10000. ) zscale=10000.
if ( zscale < 1500.  ) zscale=1500.


end  function fzscale


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

real function fx1 ( rc , renv, dz, dmfc  ) result( rc1 ) 
!===================================================
! Calculates one level of in-draft transport
!
!    ----- rc(0)
!      ++++++++++  dmfc(0)
!    ----- rc(1)              ~~  renv(1) , dz(1)
!      ++++++++++  dmfc(1)
!     
!===================================================
real, intent(inout), dimension(0:1)  :: dmfc,rc
real, intent(in) :: renv,dz             !! ,rc0
!!real  :: rc1

real :: EDc , EDf
     


IF ( DMFC(0)+DMFC(1) .ne. 0.0 ) Then

    EDC = ( DMFC(0)-DMFC(1) ) / DZ 
  
   
    IF (EDC >= 0. ) then     !  entraining layer
      rc1 = ( DMFC(0) * rc(0) - EDc*renv*DZ )/DMFC(1)
    endif     


    IF (EDC < 0. ) then     !  detraining layer
      rc1 = DMFC(0) * rc(0) / ( DMFC(1) + EDC*DZ )
    endif

else

     rc1 = rc(1)

endif

          !fx1( QVc(L-1) , QV(L-1:L), RH(L-1:L), & 
          !             DZ(L-1:L), DMFC(L-1:L) )
end  function fx1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


real function eup1 ( re , rc, dz, rho, dt, dmfc  ) result( re1 ) 
!===================================================
! Calculates one level of environmental tendency
! due to compensating mass fluxes
!
!    ----- re(0)
!      ++++++++++  dmfc(0)
!    ----- re(1)              ~~  rc(1) , dz(1), dmfxx, rho(1)
!      ++++++++++  dmfc(1)
!     
!===================================================
real, intent(inout), dimension(0:1)  :: dmfc,re
real, intent(in) ::  dz , rc , rho  , dt     !! ,rc0
!!real  :: rc1

real :: EDc , EDf, dmfxx
     

dmfxx = 0.5 * ( DMFC(0)+DMFC(1) )

EDC = -( DMFC(0)-DMFC(1) ) / DZ ! DE,DE,DEtrainment if positive 
EDC = MAX( EDC, 0. )

re1 = re(1) + (dt/rho) * ( dmfxx*( re(0)-re(1) )/Dz  + EDC*( rc - re(1) ) )
  
end  function eup1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MICRO_AA_BB(TEMP,PR,Q_SAT,AA,BB)

real, intent(in ) :: TEMP,PR,Q_SAT
real, intent(out) :: AA,BB

real  :: E_SAT

real, parameter  :: K_COND  =  2.4e-2        ! J m**-1 s**-1 K**-1
real, parameter  :: DIFFU   =  2.2e-5        ! m**2 s**-1
real, parameter  :: epsilon =  MAPL_H2OMW/MAPL_AIRMW

E_SAT = 100.* PR * Q_SAT /( EPSILON + (1.0-EPSILON)*Q_SAT )  ! (100 converts 
                                                             ! from mbar to Pa)

 AA  = ( MAPL_ALHL**2 ) / ( K_COND*MAPL_RVAP*(TEMP**2) )

 BB  = MAPL_RVAP*TEMP / ( DIFFU*(1000./PR)*E_SAT )

end subroutine MICRO_AA_BB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end MODULE DDF
