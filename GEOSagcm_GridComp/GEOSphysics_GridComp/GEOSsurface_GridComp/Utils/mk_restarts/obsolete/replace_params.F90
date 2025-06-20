PROGRAM replace_params
  implicit none

  integer :: nland_old,  nland_new
  character*400 :: tilefile
  character*400 :: old_restart, new_restart, sarithpath

  real, allocatable :: var1(:),var2(:,:)
  integer :: numrecs, allrecs, numparmrecs

  real, allocatable :: BF1(:),   BF2(:),   BF3(:),  VGWMAX(:)
  real, allocatable :: CDCR1(:), CDCR2(:), PSIS(:), BEE(:) 
  real, allocatable :: POROS(:), WPWET(:), COND(:), GNU(:)
  real, allocatable :: ARS1(:),  ARS2(:),  ARS3(:)
  real, allocatable :: ARA1(:),  ARA2(:),  ARA3(:), ARA4(:)
  real, allocatable :: ARW1(:),  ARW2(:),  ARW3(:), ARW4(:)
  real, allocatable :: TSA1(:),  TSA2(:),  TSB1(:), TSB2(:)
  real, allocatable :: ATAU2(:), BTAU2(:), ITY0(:)
  integer, allocatable :: ity_int(:)
  integer           :: ity2, nmin, type
  real, allocatable :: DP2BR(:), tmp_wgt(:,:), tmp_sum(:,:)
  real              :: zdep1, zdep2, zdep3, zmet, term1, term2
  integer           :: catindex21, catindex22, catindex23
  integer           :: catindex24, catindex25, catindex26
  real              :: frc1, frc2, rdum
  integer           :: catid, checksum, ntilesold
  integer           :: ii0, jj0, i,j,n, idum,II
  integer           :: IARGC

  
  II = iargc()

  if(II /= 4) then
     print *, "Wrong Number of arguments: ", ii
     call exit(66)
  end if

  call getarg(1,old_restart)
  call getarg(2,new_restart)
  call getarg(3,tilefile)
  call getarg(4,sarithpath)

  sarithpath   = "/land/l_data/geos5/bcs/SiB2_V2/DC/"//trim(sarithpath)

  numrecs     = 61
  numparmrecs = 30

  ! read .til file

  open (10,file=trim(tilefile),status='old',form='formatted')
  read (10,*) ntilesold
  read (10,*) 
  read (10,*) 
  read (10,*)
  read (10,*)
  read (10,*)
  read (10,*)
  read (10,*)
  nland_old=0
  do n = 1,ntilesold
     read(10,*) type
     if (type ==  100) then
        nland_old=nland_old+1
     endif
  end do
  close (10,status='keep')

  print *, ' Number of land tiles = ', nland_old

  nland_new = nland_old

  allocate (   BF1(nland_new),    BF2 (nland_new),     BF3(nland_new)  )
  allocate (VGWMAX(nland_new),   CDCR1(nland_new),   CDCR2(nland_new)  ) 
  allocate (  PSIS(nland_new),     BEE(nland_new),   POROS(nland_new)  ) 
  allocate ( WPWET(nland_new),    COND(nland_new),     GNU(nland_new)  )
  allocate (  ARS1(nland_new),    ARS2(nland_new),    ARS3(nland_new)  )
  allocate (  ARA1(nland_new),    ARA2(nland_new),    ARA3(nland_new)  )
  allocate (  ARA4(nland_new),    ARW1(nland_new),    ARW2(nland_new)  )
  allocate (  ARW3(nland_new),    ARW4(nland_new),    TSA1(nland_new)  )
  allocate (  TSA2(nland_new),    TSB1(nland_new),    TSB2(nland_new)  )
  allocate ( ATAU2(nland_new),   BTAU2(nland_new),   DP2BR(nland_new)  )
  allocate (  ITY0(nland_new), ity_int(nland_new))



  open(unit=21, file=trim(sarithpath) // '/' //'mosaic_veg_typs_fracs',form='formatted')
  open(unit=22, file=trim(sarithpath) // '/' //'bf.dat'               ,form='formatted')
  open(unit=23, file=trim(sarithpath) // '/' //'soil_param.dat'       ,form='formatted')
  open(unit=24, file=trim(sarithpath) // '/' //'ar.new'               ,form='formatted')
  open(unit=25, file=trim(sarithpath) // '/' //'ts.dat'               ,form='formatted')
  open(unit=26, file=trim(sarithpath) // '/' //'tau_param.dat'        ,form='formatted')


  print *, 'opened units'

  do n=1,nland_new 
     read (21, *) catindex21, catid, ity_int(n), ity2, frc1, frc2
     ITY0(n)=1.0*ity_int(n)

     read (22, *) catindex22, catid, GNU(n), BF1(n), BF2(n), BF3(n)

     read (23, *) catindex23, catid, idum, idum, BEE(n), PSIS(n), POROS(n), COND(n), WPWET(n), DP2BR(n)

     read (24, *) catindex24, catid, rdum, ARS1(n), ARS2(n), ARS3(n),          &
          ARA1(n), ARA2(n), ARA3(n), ARA4(n), &
          ARW1(n), ARW2(n), ARW3(n), ARW4(n)

     read (25, *) catindex25, catid, rdum, TSA1(n), TSA2(n), TSB1(n), TSB2(n)

     read (26, *) catindex26, catid, ATAU2(n), BTAU2(n), rdum, rdum

     zdep2=1000.
     zdep3=amax1(1000.,DP2BR(n))
     if (zdep2 .gt.0.75*zdep3) then
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

  close (21)
  close (22)
  close (23)
  close (24)
  close (25)
  close (26)


  print *, ' Doing restarts'

  open(unit=30, file=trim(old_restart),form='unformatted',status='old',convert='little_endian')
  open(unit=40, file=trim(new_restart),form='unformatted',status='unknown',convert='little_endian')

  allocate(var1(nland_old))
  allocate(var2(nland_old,4))

  print *, 'Opened restart files'
  print *, 30, trim(old_restart)
  print *, 40, trim(new_restart)



  write(40) BF1
  read(30) var1
  print *, "BF1",maxval(BF1), maxval(var1), minval(BF1),minval(var1)

  write(40) BF2
  read(30) var1
  print *, "BF2",maxval(BF2), maxval(var1), minval(BF2),minval(var1)

  write(40) BF3
  read(30) var1
  print *, "BF3",maxval(BF3), maxval(var1), minval(BF3),minval(var1)

  write(40) VGWMAX
  read(30) var1
  print *, "VGWMAX",maxval(VGWMAX), maxval(var1), minval(VGWMAX),minval(var1)

  write(40) CDCR1
  read(30) var1
  print *, "CDCR1",maxval(CDCR1), maxval(var1), minval(CDCR1),minval(var1)

  write(40) CDCR2
  read(30) var1
  print *, "CDCR2",maxval(CDCR2), maxval(var1), minval(CDCR2),minval(var1)

  write(40) PSIS
  read(30) var1
  print *, "PSIS",maxval(PSIS), maxval(var1), minval(PSIS),minval(var1)

  write(40) BEE
  read(30) var1
  print *, "BEE",maxval(BEE), maxval(var1), minval(BEE),minval(var1)

  write(40) POROS 
  read(30) var1
  print *, "POROS ",maxval(POROS ), maxval(var1), minval(POROS ),minval(var1)

  write(40) WPWET
  read(30) var1
  print *, "WPWET",maxval(WPWET), maxval(var1), minval(WPWET),minval(var1)

  write(40) COND
  read(30) var1
  print *, "COND",maxval(COND), maxval(var1), minval(COND),minval(var1)

  write(40) GNU
  read(30) var1
  print *, "GNU",maxval(GNU), maxval(var1), minval(GNU),minval(var1)

  write(40) ARS1
  read(30) var1
  print *, "ARS1",maxval(ARS1), maxval(var1), minval(ARS1),minval(var1)

  write(40) ARS2
  read(30) var1
  print *, "ARS2",maxval(ARS2), maxval(var1), minval(ARS2),minval(var1)

  write(40) ARS3
  read(30) var1
  print *, "ARS3",maxval(ARS3), maxval(var1), minval(ARS3),minval(var1)

  write(40) ARA1
  read(30) var1
  print *, "ARA1",maxval(ARA1), maxval(var1), minval(ARA1),minval(var1)

  write(40) ARA2
  read(30) var1
  print *, "ARA2",maxval(ARA2), maxval(var1), minval(ARA2),minval(var1)

  write(40) ARA3
  read(30) var1
  print *, "ARA3",maxval(ARA3), maxval(var1), minval(ARA3),minval(var1)

  write(40) ARA4
  read(30) var1
  print *, "ARA4",maxval(ARA4), maxval(var1), minval(ARA4),minval(var1)

  write(40) ARW1
  read(30) var1
  print *, "ARW1",maxval(ARW1), maxval(var1), minval(ARW1),minval(var1)

  write(40) ARW2
  read(30) var1
  print *, "ARW2",maxval(ARW2), maxval(var1), minval(ARW2),minval(var1)

  write(40) ARW3
  read(30) var1
  print *, "ARW3",maxval(ARW3), maxval(var1), minval(ARW3),minval(var1)

  write(40) ARW4
  read(30) var1
  print *, "ARW4",maxval(ARW4), maxval(var1), minval(ARW4),minval(var1)

  write(40) TSA1
  read(30) var1
  print *, "TSA1",maxval(TSA1), maxval(var1), minval(TSA1),minval(var1)

  write(40) TSA2
  read(30) var1
  print *, "TSA2",maxval(TSA2), maxval(var1), minval(TSA2),minval(var1)

  write(40) TSB1
  read(30) var1
  print *, "TSB1",maxval(TSB1), maxval(var1), minval(TSB1),minval(var1)

  write(40) TSB2
  read(30) var1
  print *, "TSB2",maxval(TSB2), maxval(var1), minval(TSB2),minval(var1)

  write(40) ATAU2
  read(30) var1
  print *, "ATAU2",maxval(ATAU2), maxval(var1), minval(ATAU2),minval(var1)

  write(40) BTAU2
  read(30) var1
  print *, "BTAU2",maxval(BTAU2), maxval(var1), minval(BTAU2),minval(var1)

  write(40) ITY0
  read(30) var1
  print *, "ITY0",maxval(ITY0), maxval(var1), minval(ITY0),minval(var1)


  print *, 'Wrote parameters'

  do n=1,2
     read (30) var2
     write(40) var2
  end do

  do n=1,20
     read (30) var1
     write(40) var1
  enddo

  do n=1,4
     read (30) var2
     write(40) var2
  end do

  do n=1,4
     read (30) var1
     write(40) var1
  enddo

     read (30) var2
     write(40) var2

END PROGRAM replace_params
