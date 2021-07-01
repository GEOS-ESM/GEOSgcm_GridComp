#include "MAPL_Generic.h"
PROGRAM dbg_cnlsm_offline

!
! This offline driver is good for debugging CNLSM offline.
! The first step will be re-compiling GEOS_CatchCNGridComp.F90 using 
! DBG_CNLSM_INPUTS directive and running - the model would generate 3 files:
! catchcn_inputs.data, catchcn_params.data, and catchcn_updates.data.
! The driver reads catchcn_params.data and catchcn_updates.data to initialize
! parameters and model prognostics, and continues the model integration while
! reading input variables from catchcn_inputs.data at every timestep. 
! - Sarith Mahanama (9-1-2018)
!
use ESMF
use MAPL_ConstantsMod
use MAPL_ExceptionHandling
USE SURFPARAMS
use CATCHMENT_CN_MODEL
!USE MAPL_BaseMod,      ONLY: MAPL_UNDEF
USE catch_constants,   ONLY:                &
     N_SNOW         => CATCH_N_SNOW,        &
     N_GT           => CATCH_N_GT

implicit none

type :: c_inputs
   
   real :: PCU       
   real :: PLS       
   real :: SNO  
   real :: ICE
   real :: FRZR
   real :: UUU   
   real :: EVSBTS 
   real :: DEVSBTS
   real :: DEDTCS 
   real :: SHSBTS 
   real :: DHSDQAS
   real :: DSHSBTS
   real :: EVSBTT 
   real :: DEVSBTT
   real :: DEDTCT 
   real :: SHSBTT 
   real :: DHSDQAT
   real :: DSHSBTT
   real :: EVSBTW 
   real :: DEVSBTW
   real :: DEDTCW 
   real :: SHSBTW 
   real :: DHSDQAW
   real :: DSHSBTW
   real :: EVSBTN 
   real :: DEVSBTN
   real :: DEDTCN 
   real :: SHSBTN 
   real :: DHSDQAN
   real :: DSHSBTN
   real :: TA        
   real :: QA        
   real :: RAS       
   real :: RAT       
   real :: RAW       
   real :: RASN      
   real :: ZTH       
   real :: SWNETFREE 
   real :: SWNETSNOW 
   real :: LWDNSRF   
   real :: PS        
   real :: LAI0      
   real :: GRN0      
   real :: SQSCAT    
   real :: RSL1      
   real :: RSL2      
   real :: RDC      
   real :: QSATS
   real :: DQSS 
   real :: ALW1
   real :: BLW1
   real :: QSATT
   real :: DQST 
   real :: ALW2
   real :: BLW2
   real :: QSATW
   real :: DQSW 
   real :: ALW3
   real :: BLW3
   real :: QSATN
   real :: DQSN 
   real :: ALW4
   real :: BLW4
   real :: RCSAT
   real :: DRCSDT
   real :: DRCSDQ
   real :: RCUNS
   real :: DRCUDT
   real :: DRCUDQ
   
end type c_inputs

type :: c_params
   
   real :: LONS
   real :: LATS
   real :: DZSF
   integer :: VEG1
   integer :: VEG2
   real :: FVEG1
   real :: FVEG2
   real :: BF1     
   real :: BF2     
   real :: BF3     
   real :: VGWMAX  
   real :: CDCR1   
   real :: CDCR2   
   real :: PSIS    
   real :: BEE     
   real :: POROS   
   real :: WPWET   
   real :: COND    
   real :: GNU     
   real :: ARS1    
   real :: ARS2    
   real :: ARS3    
   real :: ARA1    
   real :: ARA2    
   real :: ARA3    
   real :: ARA4    
   real :: ARW1    
   real :: ARW2    
   real :: ARW3    
   real :: ARW4    
   real :: TSA1    
   real :: TSA2    
   real :: TSB1    
   real :: TSB2    
   real :: ATAU    
   real :: BTAU     
   
end type c_params

type :: c_updates
   real :: TG1
   real :: TG2
   real :: TG4  
   real :: TC1
   real :: TC2
   real :: TC4
   real :: QC1
   real :: QC2
   real :: QC4
   real :: CAPAC
   real :: CATDEF
   real :: RZEXC
   real :: SRFEXC
   real :: GHTCNT(N_gt)
   real :: WESNN (N_snow)
   real :: HTSNNN(N_snow)
   real :: SNDZN (N_snow)
   
end type c_updates

real, parameter    :: MAPL_UNDEF = 1.0e15	
character*200 :: scratch_dir
integer :: ntiles,n,t,n_global,typ,iargc,NT
integer :: tid=0,ntsteps=10000
real :: DT,pfrac
character*20 :: cdum
character*128        :: ARG
character*1          :: opt
type (c_inputs ), dimension(:) ,allocatable :: c_input
type (c_params ), dimension(:) ,allocatable :: c_param
type (c_updates), dimension(:) ,allocatable :: c_upds

character(len=ESMF_MAXSTR) :: LAND_PARAMS ! land parameter option
INTEGER,allocatable,dimension(:) :: CAT_ID
REAL :: lonbeg,lonend,latbeg,latend
REAL, allocatable, dimension (:,:) :: WESNN  ,HTSNNN ,SNDZN ,GHTCNT
REAL, allocatable, dimension (:) ::  &
     EVAPOUT ,SHOUT   ,RUNOFF ,EVPINT ,EVPSOI ,EVPVEG ,EVPICE ,&
     BFLOW   ,RUNSURF ,SMELT  ,HLWUP  ,SWNDSRF,HLATN  ,QINFIL ,&
     AR1     ,AR2     ,RZEQ   ,GHFLX  ,TPSN   ,ASNOW          ,&
     GHFLXSNO                                                 ,&
     GHFLXTSKIN                                               ,&
     TP1     ,TP2     ,TP3    ,TP4    ,TP5    ,TP6    ,SFMC   ,&
     RZMC    ,PRMC    ,ENTOT  ,WTOT   ,WCHANGE,ECHANGE,HSNACC ,&
     EVACC   ,SHACC   ,SHSNOW                                 ,&
     AVETSNOW,WAT10CM, WATSOI ,ICESOI                         ,&
     LHSNOW1, LWUPSNOW1, LWDNSNOW1, NETSWSNOW                 ,&
     TCSORIG1, TPSN1IN1, TPSN1OUT1, LHACC, TSURF             

real, allocatable               :: TC1_0(:), TC2_0(:),  TC4_0(:)
real, allocatable               :: QA1_0(:), QA2_0(:),  QA4_0(:)

real, allocatable, dimension (:) :: DWLAND, PRLAND, EVLAND, SPWATR
real, allocatable, dimension (:) :: SUBLIM, SUMEV, EB
real, allocatable, dimension (:) :: DHLAND, SWLAND, LWLAND, SHLAND, LHLAND, SNOLAND, SPLAND, SPSNOW
real, allocatable, dimension (:,:) :: FR
logical :: is_OFFLINE
integer :: OFFLINE_MODE = 0
integer :: LSM_CHOICE = 2       ! hardcoded here as Catchment-CN (LSM_CHOICE = 2), since this routine is CN specific

lonbeg = MAPL_UNDEF
lonend = MAPL_UNDEF
latbeg = MAPL_UNDEF
latend = MAPL_UNDEF

n = iargc()

if(n < 2) then	
     print *,"Usage : ./dbg_clsm_offline  -s $SCRDIR -i TID -t ntsteps -o OFFLINE_MODE"
     print *, "-s: Scratch Directory (Mandatory)"
     print *, "-i: Tile index of the problematic tile (Optional)"
     print *, "-t: Number of timesteps (Optional)"
     call exit(1)
end if

n = 1
    call getarg(n,arg)
    do while(arg(1:1)=='-')
       opt=arg(2:2)
       if(len(trim(arg))==2) then
          if(scan(opt,'zvh')==0) then
             n = n + 1
             call getarg(n,arg)
          endif
       else
          arg = arg(3:)
       end if
       select case (opt)
       case ('s')
         scratch_dir = trim(arg)
       case ('i')
         read(arg,'(i10)')tid
       case ('t')
         read(arg,'(i10)')ntsteps
       case ('o')
         read(arg,'(i10)')OFFLINE_MODE
       end select
       n = n + 1
       call getarg(n,arg)
    end do

is_OFFLINE = OFFLINE_MODE /= 0

call ESMF_ConfigGetAttribute (SCF, label='LAND_PARAMS:', value=LAND_PARAMS, DEFAULT="CN_CLM4", __RC__ )
call SurfParams_init(LAND_PARAMS,LSM_CHOICE,rc=status) 
_VERIFY(STATUS)

open (10,file=trim(scratch_dir)//'/catchcn_inputs.data' ,form ='unformatted', &
     action ='read', status ='old',convert='little_endian')
open (11,file=trim(scratch_dir)//'/catchcn_params.data' ,form ='unformatted', & 
     action ='read', status ='old',convert='little_endian')
open (12,file=trim(scratch_dir)//'/catchcn_updates.data',form ='unformatted', &
     action ='read', status ='old',convert='little_endian')

read (11) ntiles

if(ntiles == 0) then
  
  open (15,file = trim(scratch_dir)//'/tile.data' ,form ='formatted', &
     action ='read', status ='old')
     read (15,*) n_global
     do n = 1,7
     	read (15,'(a)')cdum
     enddo

     do n = 1,n_global
       read (15,*) typ
       if(typ == 100) ntiles = n
     end do
     close (15,status='keep')
endif

read (11) dt
read (11) pfrac

print *,'======================================================'
print *,'                       '
print *,'Number of tiles : ',ntiles
print *,'HeartBeat       : ',DT
print *,'                       '
print *,'======================================================'

allocate (c_input (1:ntiles))
allocate (c_param (1:ntiles))
allocate (c_upds  (1:ntiles))
allocate (cat_id  (1:ntiles))
allocate (GHTCNT(1:N_gt  ,1:ntiles))
allocate (WESNN (1:N_snow,1:ntiles))
allocate (HTSNNN(1:N_snow,1:ntiles))
allocate (SNDZN (1:N_snow,1:ntiles))  

allocate (EVAPOUT  (1:ntiles))
allocate (SHOUT    (1:ntiles))
allocate (RUNOFF   (1:ntiles))
allocate (EVPINT   (1:ntiles))
allocate (EVPSOI   (1:ntiles))
allocate (EVPVEG   (1:ntiles))
allocate (EVPICE   (1:ntiles))
allocate (BFLOW    (1:ntiles))
allocate (RUNSURF  (1:ntiles))
allocate (SMELT    (1:ntiles))
allocate (HLWUP    (1:ntiles))
allocate (SWNDSRF  (1:ntiles))
allocate (HLATN    (1:ntiles))
allocate (QINFIL   (1:ntiles))
allocate (AR1      (1:ntiles))
allocate (AR2      (1:ntiles))
allocate (RZEQ     (1:ntiles))
allocate (GHFLX    (1:ntiles))
allocate (GHFLXSNO    (1:ntiles))
allocate (GHFLXTSKIN    (1:ntiles))
allocate (TPSN     (1:ntiles))
allocate (ASNOW    (1:ntiles))
allocate (TP1      (1:ntiles))
allocate (TP2      (1:ntiles))
allocate (TP3      (1:ntiles))
allocate (TP4      (1:ntiles))
allocate (TP5      (1:ntiles))
allocate (TP6      (1:ntiles))
allocate (SFMC     (1:ntiles))
allocate (RZMC     (1:ntiles))
allocate (PRMC     (1:ntiles))
allocate (ENTOT    (1:ntiles))
allocate (WTOT     (1:ntiles))
allocate (WCHANGE  (1:ntiles))
allocate (ECHANGE  (1:ntiles))
allocate (HSNACC   (1:ntiles))
allocate (EVACC    (1:ntiles))
allocate (SHACC    (1:ntiles))
allocate (SHSNOW   (1:ntiles))
allocate (LHACC    (1:ntiles))
allocate (TSURF    (1:ntiles))
allocate (AVETSNOW (1:ntiles))
allocate (WAT10CM  (1:ntiles))
allocate (WATSOI   (1:ntiles))
allocate (ICESOI   (1:ntiles))
allocate (LHSNOW1   (1:ntiles))
allocate (LWUPSNOW1   (1:ntiles))
allocate (LWDNSNOW1   (1:ntiles))
allocate (NETSWSNOW   (1:ntiles))
allocate (TCSORIG1   (1:ntiles))
allocate (TPSN1IN1   (1:ntiles))
allocate (TPSN1OUT1   (1:ntiles))

allocate(TC1_0    (NTILES))
allocate(TC2_0    (NTILES))
allocate(TC4_0    (NTILES))
allocate(QA1_0    (NTILES))
allocate(QA2_0    (NTILES))
allocate(QA4_0    (NTILES))

NT = NTILES
if (tid /=0) NT = 1

allocate (DWLAND   (1 : NT)) 
allocate (PRLAND   (1 : NT)) 
allocate (EVLAND   (1 : NT)) 
allocate (SPWATR   (1 : NT)) 
allocate (SUBLIM   (1 : NT)) 
allocate (DHLAND   (1 : NT)) 
allocate (SWLAND   (1 : NT)) 
allocate (LWLAND   (1 : NT)) 
allocate (SHLAND   (1 : NT)) 
allocate (LHLAND   (1 : NT)) 
allocate (SNOLAND  (1 : NT)) 
allocate (SPLAND   (1 : NT)) 
allocate (SPSNOW   (1 : NT)) 
allocate (SUMEV    (1 : NT))
allocate (EB       (1 : NT))
allocate (FR    (1 : NT, 4))

do n = 1,ntiles
   cat_id (n) = n
end do

call read_catch_params  (11,ntiles,c_param)
call read_catch_updates (12,ntiles,c_upds )

do n = 1, N_gt
    GHTCNT (n,:) = c_upds%GHTCNT(n)
end do
   
do n = 1, N_snow
      
   WESNN (n,:) = c_upds%WESNN (n)
   HTSNNN(n,:) = c_upds%HTSNNN(n)
   SNDZN (n,:) = c_upds%SNDZN (n)
      
end do


close (11,status = 'keep')
close (12,status = 'keep')

do t = 1,ntsteps

   print *, 'Processing time step : ',t

   call read_catch_inputs  (10,ntiles,c_input)
      
   if(tid == 0) then

! The model is run on the entire tile space 
! AMM - lons and lats arguments to catchment

      call  CATCHCN ( NTILES, c_param%LONS, c_param%LATS, DT, PFRAC, cat_id, &
           c_param%VEG1,c_param%VEG2,c_param%FVEG1,c_param%FVEG2,c_param%DZSF   ,&
           c_input%PCU, c_input%PLS, c_input%SNO, c_input%ICE, c_input%FRZR, c_input%UUU                   ,&
           c_input%EVSBTS ,c_input%DEVSBTS ,c_input%DEDTCS   ,&
           c_input%SHSBTS ,c_input%DHSDQAS ,c_input%DSHSBTS  ,&  
           c_input%EVSBTT ,c_input%DEVSBTT ,c_input%DEDTCT   ,&
           c_input%SHSBTT ,c_input%DHSDQAT ,c_input%DSHSBTT  ,&
           c_input%EVSBTW ,c_input%DEVSBTW ,c_input%DEDTCW   ,&
           c_input%SHSBTW ,c_input%DHSDQAW ,c_input%DSHSBTW  ,&
           c_input%EVSBTN,c_input%DEVSBTN,c_input%DEDTCN   ,  &
           c_input%SHSBTN,c_input%DHSDQAN ,c_input%DSHSBTN ,  &
           c_input%TA     ,c_input%QA      ,&
           c_input%RAS    ,c_input%RAT   ,c_input%RAW  ,c_input%RASN     ,&   
           c_input%ZTH    ,c_input%SWNETFREE, &
           c_input%SWNETSNOW ,c_input%LWDNSRF   ,&
           c_input%PS     ,c_input%LAI0  ,c_input%GRN0   ,&
           c_input%SQSCAT    ,c_input%RSL1      ,&
           c_input%RSL2   ,c_input%RDC                                     ,&
           c_input%QSATS  ,c_input%DQSS  ,c_input%ALW1  ,c_input%BLW1      ,&
           c_input%QSATT  ,c_input%DQST  ,c_input%ALW2  ,c_input%BLW2      ,&
           c_input%QSATW  ,c_input%DQSW  ,c_input%ALW3  ,c_input%BLW3      ,&
           c_input%QSATN  ,c_input%DQSN  ,c_input%ALW4  ,c_input%BLW4      ,&
           c_input%RCSAT  ,c_input%DRCSDT  ,c_input%DRCSDQ  ,c_input%RCUNS  ,c_input%DRCUDT ,c_input%DRCUDQ,&
           c_param%BF1    ,c_param%BF2   ,c_param%BF3   ,c_param%VGWMAX ,&
           c_param%CDCR1 ,c_param%CDCR2 ,c_param%PSIS	             ,&
           c_param%BEE    ,c_param%POROS ,c_param%WPWET ,c_param%COND   ,&
           c_param%GNU	                                             ,&
           c_param%ARS1   ,c_param%ARS2  ,c_param%ARS3  ,c_param%ARA1   ,&
           c_param%ARA2  ,c_param%ARA3  ,c_param%ARA4	             ,&
           c_param%ARW1   ,c_param%ARW2  ,c_param%ARW3  ,c_param%ARW4   ,&
           c_param%TSA1  ,c_param%TSA2  ,c_param%TSB1  ,c_param%TSB2   ,&
           c_param%ATAU   ,c_param%BTAU  ,.false.			     ,&
           c_upds%TG1     ,c_upds%TG2    ,c_upds%TG4		     ,& 
           c_upds%TC1     ,c_upds%TC2    ,c_upds%TC4		     ,& 
           c_upds%QC1     ,c_upds%QC2    ,c_upds%QC4		     ,&
           c_upds%CAPAC   ,c_upds%CATDEF ,c_upds%RZEXC ,c_upds%SRFEXC ,&
           GHTCNT,WESNN ,HTSNNN ,SNDZN ,&
           EVAPOUT ,SHOUT   ,RUNOFF ,EVPINT ,EVPSOI ,EVPVEG ,EVPICE ,&
           BFLOW   ,RUNSURF ,SMELT  ,HLWUP  ,SWNDSRF,HLATN  ,QINFIL ,&
           AR1     ,AR2     ,RZEQ   ,GHFLX  ,GHFLXSNO, GHFLXTSKIN, TPSN   ,ASNOW          ,&
           TP1     ,TP2     ,TP3    ,TP4    ,TP5    ,TP6    ,SFMC   ,&
           RZMC    ,PRMC    ,ENTOT  ,WTOT   ,WCHANGE,ECHANGE,HSNACC ,&
           EVACC   ,SHACC   ,TSURF, SHSNOW                                 ,&
           AVETSNOW,WAT10CM, WATSOI ,ICESOI                ,&
           LHSNOW1, LWUPSNOW1, LWDNSNOW1, NETSWSNOW             ,&
           TCSORIG1, TPSN1IN1, TPSN1OUT1,                       &
           TC1_0=TC1_0, TC2_0=TC2_0, TC4_0=TC4_0                ,&
           QA1_0=QA1_0, QA2_0=QA2_0, QA4_0=QA4_0, LHACC=LHACC )
      
!      if( .not. is_OFFLINE) then
!         !amm add correction term to latent heat diagnostics (HLATN is always allocated)
!         !    this will impact the export LHLAND
!         HLATN = HLATN - LHACC
!         ! also add some portion of the correction term to evap from soil, int, veg and snow
!         SUMEV = EVPICE+EVPSOI+EVPVEG+EVPINT
!         where(SUMEV>0.)
!            EVPICE = EVPICE - EVACC*EVPICE/SUMEV
!            EVPSOI = EVPSOI - EVACC*EVPSOI/SUMEV
!            EVPINT = EVPINT - EVACC*EVPINT/SUMEV
!            EVPVEG = EVPVEG - EVACC*EVPVEG/SUMEV
!         endwhere
!      endif
      
      FR(:,1) =           AR1  * (1-ASNOW)
      FR(:,2) =           AR2  * (1-ASNOW)
      FR(:,3) = (1.0-(AR1+AR2))* (1-ASNOW)
      FR(:,4) =                     ASNOW
      
      DWLAND = WCHANGE
      PRLAND = c_input%PCU+c_input%PLS+c_input%SNO &
            + c_input%ICE +  + c_input%FRZR
      EVLAND = EVAPOUT-EVACC
      BFLOW  = BFLOW
      SPWATR = EVACC 
      SUBLIM = EVPICE*(1./MAPL_ALHS)*FR(:,4)
      DHLAND = ECHANGE
      SWLAND = SWNDSRF 
      LWLAND = c_input%LWDNSRF - HLWUP
      SHLAND = SHOUT -SHACC 
      LHLAND = HLATN
      SNOLAND = c_input%SNO
      SPLAND = SHACC 
      SPSNOW = HSNACC 

      print *,'Mean WB : ', SUM (DWLAND - (PRLAND+ EVLAND - RUNOFF - BFLOW + SPWATR))/NTILES

      EB =  DHLAND - (SWLAND + LWLAND - SHLAND - LHLAND - MAPL_ALHF*SNOLAND - SPSNOW - SPLAND)
      print *,'Mean, Max, MAXLOC EB : ', SUM (EB) / NTILES, maxval(EB), maxloc(EB,dim=1)

   else

      ! The model is run on the problematic tile
      
      call CATCHCN ( 1, c_param(tid:tid)%lons, c_param(tid:tid)%lats                ,&
           DT	    , pfrac      , cat_id(tid:tid), c_param(tid:tid)%VEG1, c_param(tid:tid)%VEG2, &
           c_param(tid:tid)%FVEG1, c_param(tid:tid)%FVEG2, c_param(tid:tid)%DZSF    ,& 
           c_input(tid:tid)%PCU ,c_input(tid:tid)%PLS ,c_input(tid:tid)%SNO , &
           c_input(tid:tid)%ICE, c_input(tid:tid)%FRZR, c_input(tid:tid)%UUU, &
           c_input(tid:tid)%EVSBTS ,c_input(tid:tid)%DEVSBTS ,c_input(tid:tid)%DEDTCS   ,&
           c_input(tid:tid)%SHSBTS ,c_input(tid:tid)%DHSDQAS ,c_input(tid:tid)%DSHSBTS  ,&  
           c_input(tid:tid)%EVSBTT ,c_input(tid:tid)%DEVSBTT ,c_input(tid:tid)%DEDTCT   ,&
           c_input(tid:tid)%SHSBTT ,c_input(tid:tid)%DHSDQAT ,c_input(tid:tid)%DSHSBTT  ,&
           c_input(tid:tid)%EVSBTW ,c_input(tid:tid)%DEVSBTW ,c_input(tid:tid)%DEDTCW   ,&
           c_input(tid:tid)%SHSBTW ,c_input(tid:tid)%DHSDQAW ,c_input(tid:tid)%DSHSBTW  ,&
           c_input(tid:tid)%EVSBTN,c_input(tid:tid)%DEVSBTN,c_input(tid:tid)%DEDTCN   ,&
           c_input(tid:tid)%SHSBTN,c_input(tid:tid)%DHSDQAN ,c_input(tid:tid)%DSHSBTN ,&
           c_input(tid:tid)%TA     ,c_input(tid:tid)%QA      ,&
           c_input(tid:tid)%RAS    ,c_input(tid:tid)%RAT   ,c_input(tid:tid)%RAW    ,c_input(tid:tid)%RASN     ,&
           c_input(tid:tid)%ZTH    ,c_input(tid:tid)%SWNETFREE,&
           c_input(tid:tid)%SWNETSNOW ,c_input(tid:tid)%LWDNSRF   ,&
           c_input(tid:tid)%PS     ,c_input(tid:tid)%LAI0  ,c_input(tid:tid)%GRN0   ,   &
           c_input(tid:tid)%SQSCAT    ,c_input(tid:tid)%RSL1      ,&
           c_input(tid:tid)%RSL2   ,c_input(tid:tid)%RDC                                     ,&
           c_input(tid:tid)%QSATS  ,c_input(tid:tid)%DQSS  ,c_input(tid:tid)%ALW1  ,c_input(tid:tid)%BLW1      ,&
           c_input(tid:tid)%QSATT  ,c_input(tid:tid)%DQST  ,c_input(tid:tid)%ALW2  ,c_input(tid:tid)%BLW2      ,&
           c_input(tid:tid)%QSATW  ,c_input(tid:tid)%DQSW  ,c_input(tid:tid)%ALW3  ,c_input(tid:tid)%BLW3      ,&
           c_input(tid:tid)%QSATN  ,c_input(tid:tid)%DQSN  ,c_input(tid:tid)%ALW4  ,c_input(tid:tid)%BLW4      ,&
           c_input(tid:tid)%RCSAT  ,c_input(tid:tid)%DRCSDT  ,c_input(tid:tid)%DRCSDQ  ,c_input(tid:tid)%RCUNS  , &
           c_input(tid:tid)%DRCUDT ,c_input(tid:tid)%DRCUDQ     ,&
           c_param(tid:tid)%BF1    ,c_param(tid:tid)%BF2   ,c_param(tid:tid)%BF3   ,c_param(tid:tid)%VGWMAX ,&
           c_param(tid:tid)%CDCR1 ,c_param(tid:tid)%CDCR2 ,c_param(tid:tid)%PSIS	             ,&
           c_param(tid:tid)%BEE    ,c_param(tid:tid)%POROS ,c_param(tid:tid)%WPWET ,c_param(tid:tid)%COND   ,&
           c_param(tid:tid)%GNU	                                             ,&
           c_param(tid:tid)%ARS1   ,c_param(tid:tid)%ARS2  ,c_param(tid:tid)%ARS3  ,c_param(tid:tid)%ARA1   ,&
           c_param(tid:tid)%ARA2  ,c_param(tid:tid)%ARA3  ,c_param(tid:tid)%ARA4	             ,&
           c_param(tid:tid)%ARW1   ,c_param(tid:tid)%ARW2  ,c_param(tid:tid)%ARW3  ,c_param(tid:tid)%ARW4   ,&
           c_param(tid:tid)%TSA1  ,c_param(tid:tid)%TSA2  ,c_param(tid:tid)%TSB1  ,c_param(tid:tid)%TSB2   ,&
           c_param(tid:tid)%ATAU   ,c_param(tid:tid)%BTAU  ,.false.			     ,&
           c_upds(tid:tid)%TG1     ,c_upds(tid:tid)%TG2    ,c_upds(tid:tid)%TG4		     ,& 
           c_upds(tid:tid)%TC1     ,c_upds(tid:tid)%TC2    ,c_upds(tid:tid)%TC4		     ,& 
           c_upds(tid:tid)%QC1     ,c_upds(tid:tid)%QC2    ,c_upds(tid:tid)%QC4		     ,&
           c_upds(tid:tid)%CAPAC   ,c_upds(tid:tid)%CATDEF ,c_upds(tid:tid)%RZEXC ,c_upds(tid:tid)%SRFEXC ,&
           GHTCNT(1:n_gt,tid:tid)   ,WESNN(1:n_snow,tid:tid)  ,HTSNNN(1:n_snow,tid:tid) ,SNDZN(1:n_snow,tid:tid)  ,&
           EVAPOUT(tid:tid) ,SHOUT(tid:tid)   ,RUNOFF(tid:tid) ,EVPINT(tid:tid) ,EVPSOI(tid:tid) ,EVPVEG(tid:tid) ,EVPICE(tid:tid) ,&
           BFLOW(tid:tid)   ,RUNSURF(tid:tid) ,SMELT(tid:tid)  ,HLWUP(tid:tid)  ,SWNDSRF(tid:tid),HLATN(tid:tid)  ,QINFIL(tid:tid) ,&
           AR1(tid:tid)     ,AR2(tid:tid)     ,RZEQ(tid:tid)   ,GHFLX(tid:tid)  ,GHFLXSNO(tid:tid), GHFLXTSKIN(tid:tid), TPSN(tid:tid)   ,ASNOW(tid:tid)          ,&
           TP1(tid:tid)     ,TP2(tid:tid)     ,TP3(tid:tid)    ,TP4(tid:tid)    ,TP5(tid:tid)    ,TP6(tid:tid)    ,SFMC(tid:tid)   ,&
           RZMC(tid:tid)    ,PRMC(tid:tid)    ,ENTOT(tid:tid)  ,WTOT(tid:tid)   ,WCHANGE(tid:tid),ECHANGE(tid:tid),HSNACC(tid:tid) ,&
           EVACC(tid:tid)   ,SHACC(tid:tid)   ,TSURF(tid:tid) ,SHSNOW(tid:tid)                                 ,&
           AVETSNOW(tid:tid),WAT10CM(tid:tid), WATSOI(tid:tid) ,ICESOI(tid:tid)                ,&
           LHSNOW1(tid:tid), LWUPSNOW1(tid:tid), LWDNSNOW1(tid:tid), NETSWSNOW(tid:tid)             ,&
           TCSORIG1(tid:tid), TPSN1IN1(tid:tid), TPSN1OUT1(tid:tid),&
           TC1_0=TC1_0(tid:tid), TC2_0=TC2_0(tid:tid), TC4_0=TC4_0(tid:tid)                ,&
           QA1_0=QA1_0(tid:tid), QA2_0=QA2_0(tid:tid), QA4_0=QA4_0(tid:tid), LHACC=LHACC(tid:tid) )
	
!      if( .not. is_OFFLINE) then
!         !amm add correction term to latent heat diagnostics (HLATN is always allocated)
!         !    this will impact the export LHLAND
!         HLATN(tid:tid) = HLATN(tid:tid) - LHACC(tid:tid)
!         ! also add some portion of the correction term to evap from soil, int, veg and snow
!         SUMEV(:) = EVPICE(tid:tid)+EVPSOI(tid:tid)+EVPVEG+EVPINT(tid:tid)
!         if(SUMEV(1)>0.) then
!            EVPICE(tid:tid) = EVPICE(tid:tid) - EVACC(tid:tid)*EVPICE(tid:tid)/SUMEV
!            EVPSOI(tid:tid) = EVPSOI(tid:tid) - EVACC(tid:tid)*EVPSOI(tid:tid)/SUMEV
!            EVPINT(tid:tid) = EVPINT(tid:tid) - EVACC(tid:tid)*EVPINT(tid:tid)/SUMEV
!            EVPVEG(tid:tid) = EVPVEG(tid:tid) - EVACC(tid:tid)*EVPVEG(tid:tid)/SUMEV
!         endif
!      endif
      
      FR(:,1) =           AR1(tid:tid)  * (1-ASNOW(tid:tid))
      FR(:,2) =           AR2(tid:tid)  * (1-ASNOW(tid:tid))
      FR(:,3) = (1.0-(AR1(tid:tid)+AR2(tid:tid)))* (1-ASNOW(tid:tid))
      FR(:,4) =                     ASNOW(tid:tid)
      
      DWLAND(:) = WCHANGE(tid:tid)
      PRLAND(:) = c_input(tid:tid)%PCU + c_input(tid:tid)%PLS + c_input(tid:tid)%SNO &
            + c_input(tid:tid)%ICE +  + c_input(tid:tid)%FRZR
      EVLAND(:) = EVAPOUT(tid:tid)-EVACC(tid:tid)
      SPWATR(:) = EVACC(tid:tid)
      SUBLIM(:) = EVPICE(tid:tid)*(1./MAPL_ALHS)*FR(tid,4)
      DHLAND(:) = ECHANGE(tid:tid)
      SWLAND(:) = SWNDSRF(tid:tid) 
      LWLAND(:) = c_input(tid:tid)%LWDNSRF - HLWUP(tid:tid) 
      SHLAND(:) = SHOUT(tid:tid) -SHACC(tid:tid) 
      LHLAND(:) = HLATN(tid:tid) 
      SNOLAND(:) = c_input(tid:tid)%SNO
      SPLAND(:) = SHACC(tid:tid) 
      SPSNOW(:) = HSNACC(tid:tid) 
  
      print *,'WB : ', DWLAND - (PRLAND+ EVLAND - RUNOFF(tid:tid) - BFLOW(tid:tid) + SPWATR)
      print *, ' '
      print *, 'DHLAND : ',DHLAND
      print *, 'SWLAND : ',SWLAND
      print *, 'LWLAND : ',LWLAND
      print *, 'SHLAND : ',SHLAND
      print *, 'LHLAND : ',LHLAND
      print *, 'SNOLAND : ',MAPL_ALHS*SNOLAND
      print *, 'SPLAND : ',SPLAND
      print *, 'SPSNOW : ',SPSNOW

      print *,'EB : ', DHLAND - (SWLAND + LWLAND - SHLAND - LHLAND - MAPL_ALHF*SNOLAND - SPSNOW - SPLAND)
  
   endif

end do

contains

!
! -------------------------------------------------------------------
!
  subroutine read_catch_params (unit,ntiles,catchin)
    
    implicit none
    
    integer, intent (in) :: ntiles,unit
    integer :: n
    type(c_params), dimension (ntiles), intent(inout) :: catchin
    
    read (unit) catchin%LONS
    read (unit) catchin%LATS
    read (unit) catchin%DZSF
    read (unit) catchin%VEG1
    read (unit) catchin%VEG2
    read (unit) catchin%FVEG1
    read (unit) catchin%FVEG2
    read (unit) catchin%BF1     
    read (unit) catchin%BF2     
    read (unit) catchin%BF3     
    read (unit) catchin%VGWMAX  
    read (unit) catchin%CDCR1   
    read (unit) catchin%CDCR2   
    read (unit) catchin%PSIS    
    read (unit) catchin%BEE     
    read (unit) catchin%POROS   
    read (unit) catchin%WPWET   
    read (unit) catchin%COND    
    read (unit) catchin%GNU     
    read (unit) catchin%ARS1    
    read (unit) catchin%ARS2    
    read (unit) catchin%ARS3    
    read (unit) catchin%ARA1    
    read (unit) catchin%ARA2    
    read (unit) catchin%ARA3    
    read (unit) catchin%ARA4    
    read (unit) catchin%ARW1    
    read (unit) catchin%ARW2    
    read (unit) catchin%ARW3    
    read (unit) catchin%ARW4    
    read (unit) catchin%TSA1    
    read (unit) catchin%TSA2    
    read (unit) catchin%TSB1    
    read (unit) catchin%TSB2    
    read (unit) catchin%ATAU    
    read (unit) catchin%BTAU    
    
  end subroutine read_catch_params
  
  !
  ! ------------------------------------------------------
  !

  subroutine read_catch_updates (unit,ntiles,catchin)
    
    implicit none
    
    integer, intent (in) :: ntiles,unit
    type(c_updates), dimension (ntiles), intent(inout) :: catchin
    
    read (unit) catchin%TG1
    read (unit) catchin%TG2
    read (unit) catchin%TG4
    read (unit) catchin%TC1
    read (unit) catchin%TC2
    read (unit) catchin%TC4
    read (unit) catchin%QC1
    read (unit) catchin%QC2
    read (unit) catchin%QC4
    read (unit) catchin%CAPAC
    read (unit) catchin%CATDEF
    read (unit) catchin%RZEXC
    read (unit) catchin%SRFEXC
    read (unit) catchin%GHTCNT(1)
    read (unit) catchin%GHTCNT(2)
    read (unit) catchin%GHTCNT(3)
    read (unit) catchin%GHTCNT(4)
    read (unit) catchin%GHTCNT(5)
    read (unit) catchin%GHTCNT(6)
    read (unit) catchin%WESNN(1)
    read (unit) catchin%WESNN(2)
    read (unit) catchin%WESNN(3)
    read (unit) catchin%HTSNNN(1)
    read (unit) catchin%HTSNNN(2)
    read (unit) catchin%HTSNNN(3)
    read (unit) catchin%SNDZN(1)
    read (unit) catchin%SNDZN(2)
    read (unit) catchin%SNDZN(3)
    
  END subroutine read_catch_updates

!
! ------------------------------------------------------
!
  subroutine read_catch_inputs (unit,ntiles,catchin)
    
    implicit none
    
    integer, intent (in) :: ntiles,unit
    type(c_inputs), dimension (ntiles), intent(inout) :: catchin
    
    read (unit) catchin%PCU      
    read (unit) catchin%PLS      
    read (unit) catchin%SNO    
    read (unit) catchin%ICE
    read (unit) catchin%FRZR  
    read (unit) catchin%UUU      
    read (unit) catchin%EVSBTS 
    read (unit) catchin%DEVSBTS
    read (unit) catchin%DEDTCS 
    read (unit) catchin%SHSBTS 
    read (unit) catchin%DHSDQAS
    read (unit) catchin%DSHSBTS
    read (unit) catchin%EVSBTT 
    read (unit) catchin%DEVSBTT
    read (unit) catchin%DEDTCT 
    read (unit) catchin%SHSBTT 
    read (unit) catchin%DHSDQAT
    read (unit) catchin%DSHSBTT
    read (unit) catchin%EVSBTW 
    read (unit) catchin%DEVSBTW
    read (unit) catchin%DEDTCW 
    read (unit) catchin%SHSBTW 
    read (unit) catchin%DHSDQAW
    read (unit) catchin%DSHSBTW
    read (unit) catchin%EVSBTN 
    read (unit) catchin%DEVSBTN
    read (unit) catchin%DEDTCN 
    read (unit) catchin%SHSBTN 
    read (unit) catchin%DHSDQAN
    read (unit) catchin%DSHSBTN
    read (unit) catchin%TA       
    read (unit) catchin%QA       
    read (unit) catchin%RAS      
    read (unit) catchin%RAT      
    read (unit) catchin%RAW      
    read (unit) catchin%RASN     
    read (unit) catchin%ZTH      
    read (unit) catchin%SWNETFREE
    read (unit) catchin%SWNETSNOW
    read (unit) catchin%LWDNSRF  
    read (unit) catchin%PS       
    read (unit) catchin%LAI0     
    read (unit) catchin%GRN0     
    read (unit) catchin%SQSCAT   
    read (unit) catchin%RSL1     
    read (unit) catchin%RSL2     
    read (unit) catchin%RDC      
    read (unit) catchin%QSATS
    read (unit) catchin%DQSS 
    read (unit) catchin%ALW1
    read (unit) catchin%BLW1
    read (unit) catchin%QSATT
    read (unit) catchin%DQST 
    read (unit) catchin%ALW2
    read (unit) catchin%BLW2
    read (unit) catchin%QSATW
    read (unit) catchin%DQSW 
    read (unit) catchin%ALW3
    read (unit) catchin%BLW3
    read (unit) catchin%QSATN
    read (unit) catchin%DQSN 
    read (unit) catchin%ALW4
    read (unit) catchin%BLW4
    read (unit) catchin%RCSAT
    read (unit) catchin%DRCSDT
    read (unit) catchin%DRCSDQ
    read (unit) catchin%RCUNS
    read (unit) catchin%DRCUDT
    read (unit) catchin%DRCUDQ
    
  end subroutine read_catch_inputs

END PROGRAM dbg_cnlsm_offline
