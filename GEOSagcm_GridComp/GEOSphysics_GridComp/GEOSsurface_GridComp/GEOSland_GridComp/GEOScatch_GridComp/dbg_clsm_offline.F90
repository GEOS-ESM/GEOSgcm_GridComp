PROGRAM dbg_clsm_offline

!
! This offline driver is good for debugging CLSM offline.
! The first step will be re-compiling GEOS_CatchGridComp.F90 using 
! DBG_CATCH_INPUTS directive and running - the model would generate 3 files:
! catch_inputs.data, catch_params.data, and catch_updates.data.
! The driver reads catch_params.data and catch_updates.data to initialize
! parameters and model prognostics, and continues the model integration while
! reading input variables from catch_inputs.data at every timestep. 
! - Sarith Mahanama (10-10-2012)
!

use CATCHMENT_MODEL
use DBG_ROUTINES
USE MAPL_BaseMod,      ONLY: MAPL_UNDEF
USE catch_constants,   ONLY:                &
     N_SNOW         => CATCH_N_SNOW,        &
     N_GT           => CATCH_N_GT     

implicit none
character*200 :: scratch_dir
integer :: ntiles,n,t,n_global,typ,iargc
integer :: tid=0,ntsteps=10000
real :: DT,pfrac
character*20 :: cdum
real, parameter :: SFCLAY = 20.
character*128        :: ARG
character*1          :: opt
type (c_inputs ), dimension(:) ,allocatable :: c_input
type (c_params ), dimension(:) ,allocatable :: c_param
type (c_updates), dimension(:) ,allocatable :: c_upds

INTEGER,allocatable,dimension(:) :: CAT_ID
REAL, allocatable, dimension (:) ::  DZSF
REAL :: lonbeg,lonend,latbeg,latend
REAL, allocatable, dimension (:,:) :: WESNN  ,HTSNNN ,SNDZN ,GHTCNT
REAL, allocatable, dimension (:) ::  &
     EVAPOUT ,SHOUT   ,RUNOFF ,EVPINT ,EVPSOI ,EVPVEG ,EVPICE ,&
     BFLOW   ,RUNSURF ,SMELT  ,HLWUP  ,SWNDSRF,HLATN  ,QINFIL ,&
     AR1     ,AR2     ,RZEQ   ,GHFLX  ,TPSN   ,ASNOW          ,&
         GHFLXSNO                                             ,&
         GHFLXTSKIN                                           ,&
     TP1     ,TP2     ,TP3    ,TP4    ,TP5    ,TP6    ,SFMC   ,&
     RZMC    ,PRMC    ,ENTOT  ,WTOT   ,WCHANGE,ECHANGE,HSNACC ,&
     EVACC   ,SHACC   ,SHSNOW                                 ,&
     AVETSNOW,WAT10CM, WATSOI ,ICESOI                         ,&
     LHSNOW1, LWUPSNOW1, LWDNSNOW1, NETSWSNOW                 ,&
     TCSORIG1, TPSN1IN1, TPSN1OUT1                 

REAL, allocatable, dimension (:) ::  lons, lats

lonbeg = MAPL_UNDEF
lonend = MAPL_UNDEF
latbeg = MAPL_UNDEF
latend = MAPL_UNDEF

n = iargc()

if(n < 2) then	
     print *,"Usage : ./dbg_clsm_offline  -s $SCRDIR -i TID -t ntsteps"
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
       end select
       n = n + 1
       call getarg(n,arg)
    end do



open (10,file=trim(scratch_dir)//'/catch_inputs.data' ,form ='unformatted', &
     action ='read', status ='old',convert='little_endian')
open (11,file=trim(scratch_dir)//'/catch_params.data' ,form ='unformatted', & 
     action ='read', status ='old',convert='little_endian')
open (12,file=trim(scratch_dir)//'/catch_updates.data',form ='unformatted', &
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
allocate (dzsf    (1:ntiles))
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

allocate (lons  (1:ntiles))
allocate (lats  (1:ntiles))

do n = 1,ntiles
   lons (n) = 0.
   lats (n) = 0.
end do

dzsf = SFCLAY

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

   call CATCHMENT ( NTILES, lons, lats                             ,&
        DT	    , pfrac      , cat_id     , c_param%VEG, DZSF  ,& 
        c_input%PCU ,c_input%PLS ,c_input%SNO,c_input%ICE          ,&
        c_input%FRZR ,c_input%UUU                         ,&
        c_input%EVSBTS ,c_input%DEVSBTS ,c_input%TILEZERO ,&
	c_input%SHSBTS ,c_input%TILEZERO ,c_input%DSHSBTS ,&  
        c_input%EVSBTT ,c_input%DEVSBTT ,c_input%TILEZERO ,&
	c_input%SHSBTT ,c_input%TILEZERO ,c_input%DSHSBTT ,&
        c_input%EVSBTW ,c_input%DEVSBTW ,c_input%TILEZERO ,&
	c_input%SHSBTW ,c_input%TILEZERO ,c_input%DSHSBTW ,&
        c_input%EVSBTSN,c_input%DEVSBTSN,c_input%TILEZERO ,&
	c_input%SHSBTSN,c_input%TILEZERO ,c_input%DSHSBTSN,&
        c_input%TA     ,c_input%QA      ,&
        c_input%RAS    ,c_input%RAT   ,c_input%RAW    ,c_input%RASN     ,&
        c_input%ZTH    ,c_input%DRPAR ,c_input%DFPAR  ,c_input%SWNETFREE,&
	c_input%SWNETSNOW ,c_input%LWDNSRF   ,&
        c_input%PS     ,c_input%LAI0  ,c_input%GRN0   ,c_input%Z2CH     ,&
	c_input%SQSCAT    ,c_input%RSL1      ,&
        c_input%RSL2   ,c_input%RDC                                     ,&
        c_input%QSATS  ,c_input%DQSS  ,c_input%ALWX  ,c_input%BLWX      ,&
        c_input%QSATT  ,c_input%DQST  ,c_input%ALWX  ,c_input%BLWX      ,&
        c_input%QSATW  ,c_input%DQSW  ,c_input%ALWX  ,c_input%BLWX      ,&
        c_input%QSATSN ,c_input%DQSSN ,c_input%ALWX  ,c_input%BLWX      ,&
        c_param%BF1    ,c_param%BF2   ,c_param%BF3   ,c_param%VGWMAX ,&
	c_param%CDCR1 ,c_param%CDCR2 ,c_param%PSIS	             ,&
        c_param%BEE    ,c_param%POROS ,c_param%WPWET ,c_param%COND   ,&
	c_param%GNU	                                             ,&
        c_param%ARS1   ,c_param%ARS2  ,c_param%ARS3  ,c_param%ARA1   ,&
	c_param%ARA2  ,c_param%ARA3  ,c_param%ARA4	             ,&
        c_param%ARW1   ,c_param%ARW2  ,c_param%ARW3  ,c_param%ARW4   ,&
	c_param%TSA1  ,c_param%TSA2  ,c_param%TSB1  ,c_param% TSB2   ,&
        c_param%ATAU   ,c_param%BTAU  ,.false.			     ,&
        c_upds%TC1     ,c_upds%TC2    ,c_upds%TC4		     ,& 
        c_upds%QC1     ,c_upds%QC2    ,c_upds%QC4		     ,&
        c_upds%CAPAC   ,c_upds%CATDEF ,c_upds%RZEXC ,c_upds%SRFEXC ,&
	GHTCNT ,c_upds%TSURF , WESNN  ,HTSNNN ,SNDZN               ,&
        EVAPOUT ,SHOUT   ,RUNOFF ,EVPINT ,EVPSOI ,EVPVEG ,EVPICE ,&
        BFLOW   ,RUNSURF ,SMELT  ,HLWUP  ,SWNDSRF,HLATN  ,QINFIL ,&
        AR1     ,AR2     ,RZEQ   ,GHFLX  ,GHFLXSNO, GHFLXTSKIN, TPSN   ,ASNOW          ,&
        TP1     ,TP2     ,TP3    ,TP4    ,TP5    ,TP6    ,SFMC   ,&
        RZMC    ,PRMC    ,ENTOT  ,WTOT   ,WCHANGE,ECHANGE,HSNACC ,&
        EVACC   ,SHACC   ,SHSNOW                                 ,&
        AVETSNOW,WAT10CM, WATSOI ,ICESOI                         ,&
        LHSNOW1, LWUPSNOW1, LWDNSNOW1, NETSWSNOW             ,&
        TCSORIG1, TPSN1IN1, TPSN1OUT1,lonbeg,lonend,latbeg,latend                        )


   else

! The model is run on the problematic tile

   call CATCHMENT ( 1, lons(tid:tid), lats(tid:tid)                ,&
        DT	    , pfrac      , cat_id(tid:tid), c_param(tid:tid)%VEG, DZSF(tid:tid)    ,& 
        c_input(tid:tid)%PCU ,c_input(tid:tid)%PLS ,c_input(tid:tid)%SNO ,c_input(tid:tid)%UUU,&
        c_input(tid:tid)%EVSBTS ,c_input(tid:tid)%DEVSBTS ,c_input(tid:tid)%TILEZERO ,&
	c_input(tid:tid)%SHSBTS ,c_input(tid:tid)%TILEZERO ,c_input(tid:tid)%DSHSBTS ,&  
        c_input(tid:tid)%EVSBTT ,c_input(tid:tid)%DEVSBTT ,c_input(tid:tid)%TILEZERO ,&
	c_input(tid:tid)%SHSBTT ,c_input(tid:tid)%TILEZERO ,c_input(tid:tid)%DSHSBTT ,&
        c_input(tid:tid)%EVSBTW ,c_input(tid:tid)%DEVSBTW ,c_input(tid:tid)%TILEZERO ,&
	c_input(tid:tid)%SHSBTW ,c_input(tid:tid)%TILEZERO ,c_input(tid:tid)%DSHSBTW ,&
        c_input(tid:tid)%EVSBTSN,c_input(tid:tid)%DEVSBTSN,c_input(tid:tid)%TILEZERO ,&
	c_input(tid:tid)%SHSBTSN,c_input(tid:tid)%TILEZERO ,c_input(tid:tid)%DSHSBTSN,&
        c_input(tid:tid)%TA     ,c_input(tid:tid)%QA      ,&
        c_input(tid:tid)%RAS    ,c_input(tid:tid)%RAT   ,c_input(tid:tid)%RAW    ,c_input(tid:tid)%RASN     ,&
        c_input(tid:tid)%ZTH    ,c_input(tid:tid)%DRPAR ,c_input(tid:tid)%DFPAR  ,c_input(tid:tid)%SWNETFREE,&
	c_input(tid:tid)%SWNETSNOW ,c_input(tid:tid)%LWDNSRF   ,&
        c_input(tid:tid)%PS     ,c_input(tid:tid)%LAI0  ,c_input(tid:tid)%GRN0   ,c_input(tid:tid)%Z2CH     ,&
	c_input(tid:tid)%SQSCAT    ,c_input(tid:tid)%RSL1      ,&
        c_input(tid:tid)%RSL2   ,c_input(tid:tid)%RDC                                     ,&
        c_input(tid:tid)%QSATS  ,c_input(tid:tid)%DQSS  ,c_input(tid:tid)%ALWX  ,c_input(tid:tid)%BLWX      ,&
        c_input(tid:tid)%QSATT  ,c_input(tid:tid)%DQST  ,c_input(tid:tid)%ALWX  ,c_input(tid:tid)%BLWX      ,&
        c_input(tid:tid)%QSATW  ,c_input(tid:tid)%DQSW  ,c_input(tid:tid)%ALWX  ,c_input(tid:tid)%BLWX      ,&
        c_input(tid:tid)%QSATSN ,c_input(tid:tid)%DQSSN ,c_input(tid:tid)%ALWX  ,c_input(tid:tid)%BLWX      ,&
        c_param(tid:tid)%BF1    ,c_param(tid:tid)%BF2   ,c_param(tid:tid)%BF3   ,c_param(tid:tid)%VGWMAX ,&
	c_param(tid:tid)%CDCR1 ,c_param(tid:tid)%CDCR2 ,c_param(tid:tid)%PSIS	             ,&
        c_param(tid:tid)%BEE    ,c_param(tid:tid)%POROS ,c_param(tid:tid)%WPWET ,c_param(tid:tid)%COND   ,&
	c_param(tid:tid)%GNU	                                             ,&
        c_param(tid:tid)%ARS1   ,c_param(tid:tid)%ARS2  ,c_param(tid:tid)%ARS3  ,c_param(tid:tid)%ARA1   ,&
	c_param(tid:tid)%ARA2  ,c_param(tid:tid)%ARA3  ,c_param(tid:tid)%ARA4	             ,&
        c_param(tid:tid)%ARW1   ,c_param(tid:tid)%ARW2  ,c_param(tid:tid)%ARW3  ,c_param(tid:tid)%ARW4   ,&
	c_param(tid:tid)%TSA1  ,c_param(tid:tid)%TSA2  ,c_param(tid:tid)%TSB1  ,c_param(tid:tid)% TSB2   ,&
        c_param(tid:tid)%ATAU   ,c_param(tid:tid)%BTAU  ,.false.			     ,&
        c_upds(tid:tid)%TC1     ,c_upds(tid:tid)%TC2    ,c_upds(tid:tid)%TC4		     ,& 
        c_upds(tid:tid)%QC1     ,c_upds(tid:tid)%QC2    ,c_upds(tid:tid)%QC4		     ,&
        c_upds(tid:tid)%CAPAC   ,c_upds(tid:tid)%CATDEF ,c_upds(tid:tid)%RZEXC ,c_upds(tid:tid)%SRFEXC ,&
	GHTCNT(1:n_gt,tid:tid) ,c_upds(tid:tid)%TSURF , WESNN(1:n_snow,tid:tid)  ,HTSNNN(1:n_snow,tid:tid) ,SNDZN(1:n_snow,tid:tid)  ,&
        EVAPOUT(tid:tid) ,SHOUT(tid:tid)   ,RUNOFF(tid:tid) ,EVPINT(tid:tid) ,EVPSOI(tid:tid) ,EVPVEG(tid:tid) ,EVPICE(tid:tid) ,&
        BFLOW(tid:tid)   ,RUNSURF(tid:tid) ,SMELT(tid:tid)  ,HLWUP(tid:tid)  ,SWNDSRF(tid:tid),HLATN(tid:tid)  ,QINFIL(tid:tid) ,&
        AR1(tid:tid)     ,AR2(tid:tid)     ,RZEQ(tid:tid)   ,GHFLX(tid:tid)  ,GHFLXSNO(tid:tid), GHFLXTSKIN(tid:tid), TPSN(tid:tid)   ,ASNOW(tid:tid)          ,&
        TP1(tid:tid)     ,TP2(tid:tid)     ,TP3(tid:tid)    ,TP4(tid:tid)    ,TP5(tid:tid)    ,TP6(tid:tid)    ,SFMC(tid:tid)   ,&
        RZMC(tid:tid)    ,PRMC(tid:tid)    ,ENTOT(tid:tid)  ,WTOT(tid:tid)   ,WCHANGE(tid:tid),ECHANGE(tid:tid),HSNACC(tid:tid) ,&
        EVACC(tid:tid)   ,SHACC(tid:tid)   ,SHSNOW(tid:tid)                                 ,&
        AVETSNOW(tid:tid),WAT10CM(tid:tid), WATSOI(tid:tid) ,ICESOI(tid:tid)                ,&
        LHSNOW1(tid:tid), LWUPSNOW1(tid:tid), LWDNSNOW1(tid:tid), NETSWSNOW(tid:tid)             ,&
        TCSORIG1(tid:tid), TPSN1IN1(tid:tid), TPSN1OUT1(tid:tid),lonbeg,lonend,latbeg,latend                        )

endif	

   
end do

END PROGRAM dbg_clsm_offline
