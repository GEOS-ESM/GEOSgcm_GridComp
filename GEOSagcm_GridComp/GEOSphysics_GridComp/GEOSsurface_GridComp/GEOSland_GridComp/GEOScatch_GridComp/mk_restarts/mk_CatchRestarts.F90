program  mk_CatchRestarts

!  $Id$

  use MAPL_HashMod
  use MAPL_IOMod

  implicit none


  character*256 :: Usage="mk_CatchRestarts OutTileFile InTileFile InRestart SURFLAY <OutIsOld>"
  character*256 :: DataDir="OutData/clsm/"
  character*256 :: OutTileFile
  character*256 :: InTileFile
  character*256 :: InRestart
  character*256 :: DateTime
  character*256 :: OutType
  character*256 :: arg(6)

  integer :: zoom

  logical :: havedata, NewLand
  integer :: i, iargc, n, mask,j,k,ncatch,ntiles,idum,hash,last,nn,prev,Pfx
  integer :: yr,mn,yr1,mn1
  integer, pointer  :: Pf(:), Id(:), ity(:), loni(:),lono(:), lati(:), lato(:)
  integer, pointer  :: Pfi(:), Idi(:), Ido(:), idx(:), ii(:)
  real, allocatable :: lai(:), grn(:), alb(:,:)

  real, allocatable :: BF1(:),   BF2(:),   BF3(:),  VGWMAX(:)
  real, allocatable :: CDCR1(:), CDCR2(:), PSIS(:), BEE(:) 
  real, allocatable :: POROS(:), WPWET(:), COND(:), GNU(:)
  real, allocatable :: ARS1(:),  ARS2(:),  ARS3(:)
  real, allocatable :: ARA1(:),  ARA2(:),  ARA3(:), ARA4(:)
  real, allocatable :: ARW1(:),  ARW2(:),  ARW3(:), ARW4(:)
  real, allocatable :: TSA1(:),  TSA2(:),  TSB1(:), TSB2(:)
  real, allocatable :: ATAU2(:), BTAU2(:), DP2BR(:), rity(:)

  real              :: zdep1, zdep2, zdep3, zmet, term1, term2, rdum
  real, allocatable :: var1(:),var2(:,:)
  real :: timestamp(14)
  real :: SURFLAY

  logical :: InIsOld, OutIsOld

  type(MAPL_NCIO) :: InNCIO, OutNCIO
  character*256   :: vname
  character*256     :: OutFileName
  integer           :: rc
  character(len=256), parameter :: Names(29) = &
      (/'BF1   ','BF2   ','BF3   ','VGWMAX','CDCR1 ', &
        'CDCR2 ','PSIS  ','BEE   ','POROS ','WPWET ', &
        'COND  ','GNU   ','ARS1  ','ARS2  ','ARS3  ', &
        'ARA1  ','ARA2  ','ARA3  ','ARA4  ','ARW1  ', &
        'ARW2  ','ARW3  ','ARW4  ','TSA1  ','TSA2  ', & 
        'TSB1  ','TSB2  ','ATAU  ','BTAU  '/)
  logical,allocatable :: written(:)
  integer              :: ndims,filetype
  integer              :: dimSizes(3),nVars
  logical              :: file_exists
  type(MAPL_NCIO)      :: NCIOCatch  

!---------------------------------------------------------------------------

  
  I = iargc()

  if( I<5 .or. I>6 ) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     call exit(1)
  end if
 
  do n=1,I
  call getarg(n,arg(n))
  enddo
  read(arg(1),'(a)') OutTileFile
  read(arg(2),'(a)')  InTileFile
  read(arg(3),'(a)')  InRestart
  read(arg(4),*)  SURFLAY
  read(arg(5),*)  zoom

  if(I==6) then
     call getarg(6,OutType)
     OutIsOld = trim(OutType)=="OutIsOld"
  else
     OutIsOld = .false.
  endif

  if (SURFLAY.ne.20 .and. SURFLAY.ne.50) then
     print *, "You must supply a valid SURFLAY value:"
     print *, "(Ganymed-3 and earlier) SURFLAY=20.0 for Old Soil Params"
     print *, "(Ganymed-4 and later  ) SURFLAY=50.0 for New Soil Params"
     call exit(2)
  end if
  print *, 'SURFLAY: ',SURFLAY

  inquire(file=trim(DataDir)//"mosaic_veg_typs_fracs",exist=havedata)
  inquire(file=trim(DataDir)//"CLM_veg_typs_fracs"   ,exist=NewLand )

  print *, 'havedata = ',havedata

! Read Output Tile File .til file
! to get the index into the pfafsttater table

  call ReadTileFile(OutTileFile,Pf,Ido,lono,lato)

  ntiles = size(Pf)

  call ReadTileFile(InTileFile,Pfi,Idi,loni,lati)

  deallocate(Idi)
  allocate(Id(ntiles))

  call MAPL_NCIOGetFileType(InRestart, filetype,rc=rc)

  if (filetype == 0) then

     InNCIO = MAPL_NCIOOpen(InRestart,rc=rc)

     call MAPL_NCIOChangeRes(InNCIO,OutNCIO,tileSize=ntiles,rc=rc)
     i = index(InRestart,'/',back=.true.)
     OutFileName = "OutData/"//trim(InRestart(i+1:))
     call MAPL_NCIOSet( OutNCIO,filename=OutFileName )
     call MAPL_NCIOCreateFile(OutNCIO)
     call MAPL_NCIOGetDimSizes(InNCIO,nVars=nVars)
     allocate(written(nvars))
     written=.false.

  else

     open(unit=50,FILE=InRestart,form='unformatted',&
          status='old',convert='little_endian')

     do i=1,58
       read(50,end=2001)
     end do
 2001 continue
     InIsOld = I==59

     rewind(50)

     i = index(InRestart,'/',back=.true.)

     open(unit=40,FILE="OutData/"//trim(InRestart(i+1:)),form='unformatted',&
          status='unknown',convert='little_endian')

  end if

  call GetIds(loni,lati,lono,lato,Id)

  HAVE: if(havedata) then

     print *,'Working from Sariths data pretiled for this resolution'

! Get number of catchments

  open(unit=22, &
       file=trim(DataDir)//"catchment.def",status='old',form='formatted')

  read (22, *) ncatch

  close(22)

  if(ncatch==size(ido)) then
     print *, "Read ",Ncatch," land tiles."
     do i=1,ncatch
        ido(i) = i
     enddo
  else
     print *, "Number of tiles in data, ",Ncatch," does not match number in til file ", size(Ido)
     call exit(1)
  endif

  allocate(ity(ncatch),rity(ncatch))
  allocate (   BF1(ncatch),    BF2 (ncatch),     BF3(ncatch)  )
  allocate (VGWMAX(ncatch),   CDCR1(ncatch),   CDCR2(ncatch)  ) 
  allocate (  PSIS(ncatch),     BEE(ncatch),   POROS(ncatch)  ) 
  allocate ( WPWET(ncatch),    COND(ncatch),     GNU(ncatch)  )
  allocate (  ARS1(ncatch),    ARS2(ncatch),    ARS3(ncatch)  )
  allocate (  ARA1(ncatch),    ARA2(ncatch),    ARA3(ncatch)  )
  allocate (  ARA4(ncatch),    ARW1(ncatch),    ARW2(ncatch)  )
  allocate (  ARW3(ncatch),    ARW4(ncatch),    TSA1(ncatch)  )
  allocate (  TSA2(ncatch),    TSB1(ncatch),    TSB2(ncatch)  )
  allocate ( ATAU2(ncatch),    BTAU2(ncatch),   DP2BR(ncatch)  )

  inquire(file = trim(DataDir)//'/catch_params.nc4', exist=file_exists)

  if(file_exists) then
         print *,'FILE FORMAT FOR LAND BCS IS NC4'
       NCIOCatch   = MAPL_NCIOOpen(trim(DataDir)//'/catch_params.nc4',rc=rc) 
       call MAPL_VarRead ( NCIOCatch ,'OLD_ITY', rity)
       call MAPL_VarRead ( NCIOCatch ,'ARA1', ARA1)
       call MAPL_VarRead ( NCIOCatch ,'ARA2', ARA2)
       call MAPL_VarRead ( NCIOCatch ,'ARA3', ARA3)
       call MAPL_VarRead ( NCIOCatch ,'ARA4', ARA4)
       call MAPL_VarRead ( NCIOCatch ,'ARS1', ARS1)
       call MAPL_VarRead ( NCIOCatch ,'ARS2', ARS2)
       call MAPL_VarRead ( NCIOCatch ,'ARS3', ARS3)
       call MAPL_VarRead ( NCIOCatch ,'ARW1', ARW1)
       call MAPL_VarRead ( NCIOCatch ,'ARW2', ARW2)
       call MAPL_VarRead ( NCIOCatch ,'ARW3', ARW3)
       call MAPL_VarRead ( NCIOCatch ,'ARW4', ARW4)

       if( SURFLAY.eq.20.0 ) then
          call MAPL_VarRead ( NCIOCatch ,'ATAU2', ATAU2)
          call MAPL_VarRead ( NCIOCatch ,'BTAU2', BTAU2)
       endif

       if( SURFLAY.eq.50.0 ) then
          call MAPL_VarRead ( NCIOCatch ,'ATAU5', ATAU2)
          call MAPL_VarRead ( NCIOCatch ,'BTAU5', BTAU2)
       endif

       call MAPL_VarRead ( NCIOCatch ,'PSIS', PSIS)
       call MAPL_VarRead ( NCIOCatch ,'BEE', BEE)
       call MAPL_VarRead ( NCIOCatch ,'BF1', BF1)
       call MAPL_VarRead ( NCIOCatch ,'BF2', BF2)
       call MAPL_VarRead ( NCIOCatch ,'BF3', BF3)
       call MAPL_VarRead ( NCIOCatch ,'TSA1', TSA1)
       call MAPL_VarRead ( NCIOCatch ,'TSA2', TSA2)
       call MAPL_VarRead ( NCIOCatch ,'TSB1', TSB1)
       call MAPL_VarRead ( NCIOCatch ,'TSB2', TSB2)
       call MAPL_VarRead ( NCIOCatch ,'COND', COND)
       call MAPL_VarRead ( NCIOCatch ,'GNU', GNU)
       call MAPL_VarRead ( NCIOCatch ,'WPWET', WPWET)
       call MAPL_VarRead ( NCIOCatch ,'DP2BR', DP2BR)
       call MAPL_VarRead ( NCIOCatch ,'POROS', POROS)
       call MAPL_NCIOClose (NCIOCatch  )   
  else
     open(unit=21, file=trim(DataDir)//"mosaic_veg_typs_fracs",status='old',form='formatted')
     open(unit=22, file=trim(DataDir)//'bf.dat'               ,form='formatted')
     open(unit=23, file=trim(DataDir)//'soil_param.dat'       ,form='formatted')
     open(unit=24, file=trim(DataDir)//'ar.new'               ,form='formatted')
     open(unit=25, file=trim(DataDir)//'ts.dat'               ,form='formatted')
     open(unit=26, file=trim(DataDir)//'tau_param.dat'        ,form='formatted')
     
     do n=1,ncatch
        read (21,*) I, j, ITY(N)
        read (22, *) i,j, GNU(n), BF1(n), BF2(n), BF3(n)
        
        read (23, *) i,j, idum, idum, BEE(n), PSIS(n),&
             POROS(n), COND(n), WPWET(n), DP2BR(n)
        
        read (24, *) i,j, rdum, ARS1(n), ARS2(n), ARS3(n),          &
             ARA1(n), ARA2(n), ARA3(n), ARA4(n), &
             ARW1(n), ARW2(n), ARW3(n), ARW4(n)
        
        read (25, *) i,j, rdum, TSA1(n), TSA2(n), TSB1(n), TSB2(n)
        
        if( SURFLAY.eq.20.0 ) read (26, *) i,j, ATAU2(n), BTAU2(n), rdum, rdum   ! for old soil params
        if( SURFLAY.eq.50.0 ) read (26, *) i,j, rdum , rdum, ATAU2(n), BTAU2(n)  ! for new soil params
     end do

     rity = float(ity)
     CLOSE (21, STATUS = 'KEEP')
     CLOSE (22, STATUS = 'KEEP')
     CLOSE (23, STATUS = 'KEEP')
     CLOSE (24, STATUS = 'KEEP')
     CLOSE (25, STATUS = 'KEEP')

  endif

  do n=1,ncatch
    
     zdep2=1000.
     zdep3=amax1(1000.,DP2BR(n))

     if (zdep2 > 0.75*zdep3) then
        zdep2  =  0.75*zdep3              
     end if

     zdep1=20.
     zmet=zdep3/1000.

     term1=-1.+((PSIS(n)-zmet)/PSIS(n))**((BEE(n)-1.)/BEE(n))
     term2=PSIS(n)*BEE(n)/(BEE(n)-1)

     VGWMAX(n) = POROS(n)*zdep2   
     CDCR1(n)  = 1000.*POROS(n)*(zmet-(-term2*term1))   
     CDCR2(n)  = (1.-WPWET(n))*POROS(n)*zdep3
  enddo


 if (filetype /=0) then
        do i=1,30
           read(50)
        enddo
     end if

     idx => ido

else

     print *,'Working from restarts alone'

   ncatch = size(pfi)

  allocate(   rity(ncatch))
  allocate(   BF1(ncatch),    BF2 (ncatch),     BF3(ncatch)  )
  allocate (VGWMAX(ncatch),   CDCR1(ncatch),   CDCR2(ncatch)  ) 
  allocate (  PSIS(ncatch),     BEE(ncatch),   POROS(ncatch)  ) 
  allocate ( WPWET(ncatch),    COND(ncatch),     GNU(ncatch)  )
  allocate (  ARS1(ncatch),    ARS2(ncatch),    ARS3(ncatch)  )
  allocate (  ARA1(ncatch),    ARA2(ncatch),    ARA3(ncatch)  )
  allocate (  ARA4(ncatch),    ARW1(ncatch),    ARW2(ncatch)  )
  allocate (  ARW3(ncatch),    ARW4(ncatch),    TSA1(ncatch)  )
  allocate (  TSA2(ncatch),    TSB1(ncatch),    TSB2(ncatch)  )
  allocate ( ATAU2(ncatch),    BTAU2(ncatch),   DP2BR(ncatch)  )

  if (filetype == 0) then

     call MAPL_VarRead(InNCIO,names(1),BF1)
     call MAPL_VarRead(InNCIO,names(2),BF2)
     call MAPL_VarRead(InNCIO,names(3),BF3)
     call MAPL_VarRead(InNCIO,names(4),VGWMAX)
     call MAPL_VarRead(InNCIO,names(5),CDCR1)
     call MAPL_VarRead(InNCIO,names(6),CDCR2)
     call MAPL_VarRead(InNCIO,names(7),PSIS)
     call MAPL_VarRead(InNCIO,names(8),BEE)
     call MAPL_VarRead(InNCIO,names(9),POROS)
     call MAPL_VarRead(InNCIO,names(10),WPWET)

     call MAPL_VarRead(InNCIO,names(11),COND)
     call MAPL_VarRead(InNCIO,names(12),GNU)
     call MAPL_VarRead(InNCIO,names(13),ARS1)
     call MAPL_VarRead(InNCIO,names(14),ARS2)
     call MAPL_VarRead(InNCIO,names(15),ARS3)
     call MAPL_VarRead(InNCIO,names(16),ARA1)
     call MAPL_VarRead(InNCIO,names(17),ARA2)
     call MAPL_VarRead(InNCIO,names(18),ARA3)
     call MAPL_VarRead(InNCIO,names(19),ARA4)
     call MAPL_VarRead(InNCIO,names(20),ARW1)

     call MAPL_VarRead(InNCIO,names(21),ARW2)
     call MAPL_VarRead(InNCIO,names(22),ARW3)
     call MAPL_VarRead(InNCIO,names(23),ARW4)
     call MAPL_VarRead(InNCIO,names(24),TSA1)
     call MAPL_VarRead(InNCIO,names(25),TSA2)
     call MAPL_VarRead(InNCIO,names(26),TSB1)
     call MAPL_VarRead(InNCIO,names(27),TSB2)
     call MAPL_VarRead(InNCIO,names(28),ATAU2)
     call MAPL_VarRead(InNCIO,names(29),BTAU2)
     call MAPL_VarRead(InNCIO,'OLD_ITY',rITY)
 
  else

     read(50) BF1
     read(50) BF2
     read(50) BF3
     read(50) VGWMAX
     read(50) CDCR1
     read(50) CDCR2
     read(50) PSIS
     read(50) BEE
     read(50) POROS
     read(50) WPWET

     read(50) COND
     read(50) GNU
     read(50) ARS1
     read(50) ARS2
     read(50) ARS3
     read(50) ARA1
     read(50) ARA2
     read(50) ARA3
     read(50) ARA4
     read(50) ARW1

     read(50) ARW2
     read(50) ARW3
     read(50) ARW4
     read(50) TSA1
     read(50) TSA2
     read(50) TSB1
     read(50) TSB2
     read(50) ATAU2
     read(50) BTAU2
     read(50) rITY

  end if

  idx => id

  endif HAVE

  if (filetype == 0) then

     call MAPL_VarWrite(OutNCIO,names(1),BF1(Idx))
     call MAPL_VarWrite(OutNCIO,names(2),BF2(Idx))
     call MAPL_VarWrite(OutNCIO,names(3),BF3(Idx))
     call MAPL_VarWrite(OutNCIO,names(4),VGWMAX(Idx))
     call MAPL_VarWrite(OutNCIO,names(5),CDCR1(Idx))
     call MAPL_VarWrite(OutNCIO,names(6),CDCR2(Idx))
     call MAPL_VarWrite(OutNCIO,names(7),PSIS(Idx))
     call MAPL_VarWrite(OutNCIO,names(8),BEE(Idx))
     call MAPL_VarWrite(OutNCIO,names(9),POROS(Idx))
     call MAPL_VarWrite(OutNCIO,names(10),WPWET(Idx))
     call MAPL_VarWrite(OutNCIO,names(11),COND(Idx))
     call MAPL_VarWrite(OutNCIO,names(12),GNU(Idx))
     call MAPL_VarWrite(OutNCIO,names(13),ARS1(Idx))
     call MAPL_VarWrite(OutNCIO,names(14),ARS2(Idx))
     call MAPL_VarWrite(OutNCIO,names(15),ARS3(Idx))
     call MAPL_VarWrite(OutNCIO,names(16),ARA1(Idx))
     call MAPL_VarWrite(OutNCIO,names(17),ARA2(Idx))
     call MAPL_VarWrite(OutNCIO,names(18),ARA3(Idx))
     call MAPL_VarWrite(OutNCIO,names(19),ARA4(Idx))
     call MAPL_VarWrite(OutNCIO,names(20),ARW1(Idx))
     call MAPL_VarWrite(OutNCIO,names(21),ARW2(Idx))
     call MAPL_VarWrite(OutNCIO,names(22),ARW3(Idx))
     call MAPL_VarWrite(OutNCIO,names(23),ARW4(Idx))
     call MAPL_VarWrite(OutNCIO,names(24),TSA1(Idx))
     call MAPL_VarWrite(OutNCIO,names(25),TSA2(Idx))
     call MAPL_VarWrite(OutNCIO,names(26),TSB1(Idx))
     call MAPL_VarWrite(OutNCIO,names(27),TSB2(Idx))
     call MAPL_VarWrite(OutNCIO,names(28),ATAU2(Idx))
     call MAPL_VarWrite(OutNCIO,names(29),BTAU2(Idx))
     call MAPL_VarWrite(OutNCIO,'OLD_ITY',rity(Idx))

     call MAPL_NCIOGetDimSizes(InNCIO,nVars=nVars)
     do i=1,nvars
        call MAPL_NCIOGetVarName(InNCIO,i,vname)
        do j=1,29
           if ( trim(vname) == trim(names(j)) ) written(i) = .true.
        enddo
        if (trim(vname) == "OLD_ITY" ) written(i) = .true.
     enddo

     allocate(var1(size(pfi)))

     do n=1,nVars

        if (.not.written(n) ) then

           call MAPL_NCIOGetVarName(InNCIO,n,vname)

           call MAPL_NCIOVarGetDims(InNCIO,vname,nDims,dimSizes)
           if (ndims == 1) then
              call MAPL_VarRead(InNCIO,vname,var1)
              call MAPL_VarWrite(OutNCIO,vname,var1(id))
           else if (ndims == 2) then

              do j=1,dimSizes(2)
                 call MAPL_VarRead(InNCIO,vname,var1,offset1=j)
                 call MAPL_VarWrite(OutNCIO,vname,var1(id),offset1=j)
              enddo
           else if (ndims == 3) then

              do k=1,dimSizes(3)
                 do j=1,dimSizes(2)
                    call MAPL_VarRead(InNCIO,vname,var1,offset1=j,offset2=k)
                    call MAPL_VarWrite(OutNCIO,vname,var1(id),offset1=j,offset2=k)
                 enddo
              enddo

           end if

        end if

     enddo

  else

     write(40) BF1(Idx)
     write(40) BF2(Idx)
     write(40) BF3(Idx)
     write(40) VGWMAX(Idx)
     write(40) CDCR1(Idx)
     write(40) CDCR2(Idx)
     write(40) PSIS(Idx)
     write(40) BEE(Idx)
     write(40) POROS (Idx)
     write(40) WPWET(Idx)
     write(40) COND(Idx)
     write(40) GNU(Idx)
     write(40) ARS1(Idx)
     write(40) ARS2(Idx)
     write(40) ARS3(Idx)
     write(40) ARA1(Idx)
     write(40) ARA2(Idx)
     write(40) ARA3(Idx)
     write(40) ARA4(Idx)
     write(40) ARW1(Idx)
     write(40) ARW2(Idx)
     write(40) ARW3(Idx)
     write(40) ARW4(Idx)
     write(40) TSA1(Idx)
     write(40) TSA2(Idx)
     write(40) TSB1(Idx)
     write(40) TSB2(Idx)
     write(40) ATAU2(Idx)
     write(40) BTAU2(Idx)
     write(40) rITY(Idx)


     allocate(var1(size(pfi)))
     allocate(var2(size(pfi),4))

   ! TC QC

     do n=1,2
        read (50) var2
        write(40) ((var2(id(i),j),i=1,ntiles),j=1,4)
     end do

   !CAPAC CATDEF RZEXC SRFEXC ... SNDZN3

     do n=1,20
        read (50) var1
        write(40) var1(Id)
     enddo

   ! CH CM CQ FR

     do n=1,4
        read (50) var2
        write(40) ((var2(id(i),j),i=1,ntiles),j=1,4)
     end do

   ! These are the 2 prev/next pairs that dont are not
   ! in the internal in fortuna-2_0 and later. Earlier the 
   ! record are there, but their values are not needed, since
   ! they are initialized on start-up.

     if(InIsOld) then
        do n=1,4
           read (50)
        enddo
     endif

     if(OutIsOld) then
        var1 = 0.0
        do n=1,4
           write(40) var1(1:ntiles)
        end do
     endif

   ! WW

     read (50) var2
     write(40) ((var2(id(i),j),i=1,ntiles),j=1,4)

  end if

contains

#include "getids.H"

end program mk_CatchRestarts

