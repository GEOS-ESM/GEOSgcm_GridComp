PROGRAM mk_catch_internal
implicit none

integer :: im_gcm_old, jm_gcm_old
integer :: im_ocn_old, jm_ocn_old
integer :: im_gcm_new, jm_gcm_new
integer :: im_ocn_new, jm_ocn_new
integer :: ntiles_old, ntiles_new
integer :: nland_old,  nland_new

integer     qtile
parameter ( qtile = 45848 )

real,    allocatable :: lats_old(:), lats_new(:), lats_tmp(:)
real,    allocatable :: lons_old(:), lons_new(:), lons_tmp(:)
integer, allocatable ::   ii_old(:),   ii_new(:),   ii_tmp(:)
integer, allocatable ::   jj_old(:),   jj_new(:),   jj_tmp(:)
real,    allocatable ::   fr_old(:),   fr_new(:),   fr_tmp(:)
integer, allocatable ::  typ_tmp(:)

character*20  :: version1, version2
character*400 :: landir,wrkdir, old_tilefile, new_tilefile, arch, flag
character*400 :: old_rslv, old_dateline, oldtilnam
character*400 :: new_rslv, new_dateline, newtildir, newtilnam
character*400 :: old_restart, new_restart, sarithpath, home
character*400 :: maxtilnam, maxtildir, logfile
character*400 :: old_diag_grids, new_diag_grids

logical :: maxoldtoggle, twotiles
integer :: ierr, indr1, indr2, indr3, ig, jg, indx_dum, ip1, ip2
real :: fr_ocn, rdum
integer :: dum,n,nn,nta,v,loc,idum
character*4 :: bak=char(8)//char(8)//char(8)//char(8)
real, allocatable :: oldprogvars(:,:), oldparmvars(:,:), oldallvars(:,:)
real, allocatable :: newallvars(:,:)
real, allocatable :: oldvargrids(:,:,:)
real, allocatable :: newvargrids(:,:,:), dumtile(:,:), dumgrid(:,:)
real, allocatable :: tiletilevar(:,:), tilevar(:)
integer :: numrecs, allrecs, numparmrecs, numsubtiles
logical, allocatable :: ttlookup(:)
integer, allocatable :: corners_lookup(:,:)
real, allocatable :: weights_lookup(:,:)
real, allocatable :: BF1(:),   BF2(:),   BF3(:),  VGWMAX(:)
real, allocatable :: CDCR1(:), CDCR2(:), PSIS(:), BEE(:) 
real, allocatable :: POROS(:), WPWET(:), COND(:), GNU(:)
real, allocatable :: ARS1(:),  ARS2(:),  ARS3(:)
real, allocatable :: ARA1(:),  ARA2(:),  ARA3(:), ARA4(:)
real, allocatable :: ARW1(:),  ARW2(:),  ARW3(:), ARW4(:)
real, allocatable :: TSA1(:),  TSA2(:),  TSB1(:), TSB2(:)
real, allocatable :: ATAU2(:), BTAU2(:), ITY0(:)
integer, allocatable :: ity_int(:)
integer           :: ity2, nmin
real             :: frc1, frc2, dmin, dist
real, allocatable :: DP2BR(:), tmp_wgt(:,:), tmp_sum(:,:)
real              :: zdep1, zdep2, zdep3, zmet, term1, term2
integer           :: catindex21, catindex22, catindex23
integer           :: catindex24, catindex25, catindex26
integer           :: catid, checksum
integer           :: ii0, jj0, i,j
real              :: fr0, val0, ESMF_MISSING
real              :: lata, latb,lona, lonb, vaa, vbb, vab, vba
real              :: lat00_old, lon00_old, dx_old, dy_old
real              :: lonIM_old, lon0, lat0, d1,d2,d3,d4
integer           :: ia, ib, ja, jb
real              :: waa, wab, wba, wbb, wsum, tol, tempval
!real              :: mindist, olddist, thislat, thislon

! -------------------------------------------------------------------------------
! Strategy: 
! 1. Read in the "old" .til definition file
! 2. Read in the "new" .til definition file
! 3. Read in the "old" restart from a previous run
! 4. Read in Sarith's tilespace catchment parameters
! 5. Convert the prognostic variables from the old restart
!    to the new catchment definitions
!    a. Create aggregate imxjm grid of progs from old restart 
!    b. Create reasonable interolated values based on centroids
!       in the new .til definitions file.
! 6. Write the restart using stored values
! -------------------------------------------------------------------------------

! user parameters
! ---------------

  call getenv ('ARCH'        ,arch        )
  call getenv ('LANDIR'      ,landir      )
  call getenv ('WRKDIR'      ,wrkdir      )

  call getenv ('old_rslv'    ,old_rslv    )
  call getenv ('old_dateline',old_dateline)
  call getenv ('old_tilefile',old_tilefile)
  call getenv ('old_restart' ,old_restart )

  call getenv ('new_rslv'    ,new_rslv    )
  call getenv ('new_dateline',new_dateline)
  call getenv ('new_tilefile',new_tilefile)

  if( ARCH == 'OSF1'   ) flag = 'no'
  if( ARCH == 'IRIX64' ) flag = 'yes'

old_restart  = trim(wrkdir) // '/' // trim(old_restart)
new_restart  = trim(old_restart) // '.' // trim(new_rslv) // '_' // trim(new_dateline)
 sarithpath  = trim(landir) // '/'

numsubtiles =  4
numrecs     = 61  ! number of records in the restart (includes tiletile vars)
numparmrecs = 30  ! number of parameters at the beginning of restart 

allrecs = numrecs + 7*3  ! all records, with tile-tile prognostic variables expanded

allocate(ttlookup(numrecs))
ttlookup(:)     = .false. ! ttlookup specifies which records in restart are tile-tile
ttlookup(31:32) = .true.  ! or tileonly (.false.=tileonly)
ttlookup(53:56) = .true.  
ttlookup(61)    = .true.

twotiles=.false.        ! set this to true to force a two tile test
ESMF_MISSING=-999.0     ! missing value in old restart prognostic variables
tol=1.0E-6
logfile='mk_catch_restart.log'

old_diag_grids = 'old_grids.dat'
new_diag_grids = 'new_grids.dat'

! -------------------------------------------------------------------------------
! 1. Read in the old .til file and store the I, J, FR's
! -------------------------------------------------------------------------------

newtildir  = trim(sarithpath) // trim(new_dateline) // '/FV_' // trim(new_rslv)

oldtilnam  = trim(wrkdir) // '/' // trim(old_tilefile)
newtilnam  = trim(wrkdir) // '/' // trim(new_tilefile)
print *, 'newtilenam1 = ',newtilnam
!newtilnam  = trim(newtildir) // '/FV_' // trim(new_rslv) //'_'//trim(new_dateline)//'_360x180_DE_NO_TINY.til'
!newtilnam  = trim(newtildir) // '/FV_' // trim(new_rslv) //'_'//trim(new_dateline)//'_576x540_DE_NO_TINY.til'
print *, 'newtilenam2 = ',newtilnam

open(9, file=trim(logfile),action='write',form='formatted')

write (*,*)
write (*,*) '---------------------------------------------------------------------'
write (*,*) 'Reading source (old) tile definitions from:'
write (9,*) '---------------------------------------------------------------------'
write (9,*) 'Reading source (old) tile definitions from:'
write (9,*) trim(oldtilnam)
write (*,*) trim(oldtilnam)

open (10,file=trim(oldtilnam),status='old', action='read',form='formatted')
read (10,*) ntiles_old
read (10,*) dum
read (10,'(a)')version1
read (10,*)im_gcm_old
read (10,*)jm_gcm_old
read (10,'(a)')version2
read (10,*) im_ocn_old
read (10,*) jm_ocn_old
write(9,*) 'Header: ', ntiles_old, dum, trim(version1), im_gcm_old, jm_gcm_old, &
          trim(version2), im_ocn_old, jm_ocn_old

allocate(lats_tmp(ntiles_old))
allocate(lons_tmp(ntiles_old))
allocate(  fr_tmp(ntiles_old))
allocate(  ii_tmp(ntiles_old))
allocate(  jj_tmp(ntiles_old))
allocate( typ_tmp(ntiles_old))

write(*, 40, advance=trim(flag))  
nland_old=0
do n = 1,ntiles_old
   read(10,'(i10,i9,2f10.4,2i5,f10.6,3i8,f10.6,i8)',IOSTAT=ierr)typ_tmp(n),&
       indr1,lons_tmp(n),lats_tmp(n),ii_tmp(n),jj_tmp(n),fr_tmp(n),indx_dum,indr2,dum,fr_ocn,indr3
   if (typ_tmp(n) ==  100) then
      ip2=n
      nland_old=nland_old+1
   endif
   if (typ_tmp(n) == 0) then 
      ip1=n
   endif
   if(ierr /= 0) write (*,*) 'Problem reading'
   write(*, 50, advance=trim(flag)) bak, floor(float(n)/float(ntiles_old)*100)
end do
close (10,status='keep')

write(9,*) 'Last ocean index:', ip1
write(9,*) 'Last land index:', ip2
write(9,*) 'NTILES LAND:', nland_old
write(*,*)

!write(*,*) 'Packing land coordinate arrays...'

allocate ( lats_old(nland_old) )
allocate ( lons_old(nland_old) )
allocate (   fr_old(nland_old) )
allocate (   ii_old(nland_old) )
allocate (   jj_old(nland_old) )

lats_old = pack(lats_tmp, mask=typ_tmp .eq. 100)
lons_old = pack(lons_tmp, mask=typ_tmp .eq. 100)
  fr_old = pack(  fr_tmp, mask=typ_tmp .eq. 100)
  ii_old = pack(  ii_tmp, mask=typ_tmp .eq. 100)
  jj_old = pack(  jj_tmp, mask=typ_tmp .eq. 100)

write(9,*) 'lats',  size(lats_old), minval(lats_old), maxval(lats_old)
write(9,*) 'lons',  size(lons_old), minval(lons_old), maxval(lons_old)
write(9,*) 'fr  ',  size  (fr_old), minval  (fr_old), maxval  (fr_old) 
write(9,*) 'ii  ',  size  (ii_old), minval  (ii_old), maxval  (ii_old)
write(9,*) 'jj  ',  size  (jj_old), minval  (jj_old), maxval  (jj_old)

deallocate(lats_tmp)
deallocate(lons_tmp)
deallocate(  fr_tmp)
deallocate(  ii_tmp)
deallocate(  jj_tmp)
deallocate( typ_tmp)

! -------------------------------------------------------------------------------
! 2. Read in the new .til file and store the I, J, FR's
! -------------------------------------------------------------------------------

write (*,*)
write (*,*) '---------------------------------------------------------------------'
write (*,*) 'Reading source (new) tile definitions from:'
write (*,*) trim(newtilnam)
write (9,*) '---------------------------------------------------------------------'
write (9,*) 'Reading source (new) tile definitions from:'
write (9,*) trim(newtilnam)

open (10,file=trim(newtilnam),status='old',action='read',form='formatted')
read (10,*) ntiles_new
read (10,*) dum
read (10,'(a)')version1
read (10,*)im_gcm_new
read (10,*)jm_gcm_new
read (10,'(a)')version2
read (10,*) im_ocn_new
read (10,*) jm_ocn_new
write(9,*) 'Header: ', ntiles_new, dum, trim(version1), im_gcm_new, jm_gcm_new, &
          trim(version2), im_ocn_new, jm_ocn_new

allocate ( lats_tmp(ntiles_new) )
allocate ( lons_tmp(ntiles_new) )
allocate (   fr_tmp(ntiles_new) )
allocate (   ii_tmp(ntiles_new) )
allocate (   jj_tmp(ntiles_new) )
allocate (  typ_tmp(ntiles_new) )

 write(*, 40, advance=trim(flag))  
nland_new=0
do n = 1,ntiles_new
   read(10,'(i10,i9,2f10.4,2i5,f10.6,3i8,f10.6,i8)',IOSTAT=ierr)typ_tmp(n),&
       indr1,lons_tmp(n),lats_tmp(n),ii_tmp(n),jj_tmp(n),fr_tmp(n),indx_dum,indr2,dum,fr_ocn,indr3
   if (typ_tmp(n) ==  100) then
      ip2=n
      nland_new=nland_new+1
   endif
   if (typ_tmp(n) == 0) then 
      ip1=n
   endif
   if(ierr /= 0) write (*,*) 'Problem reading'
   write(*, 50, advance=trim(flag)) bak, floor(float(n)/float(ntiles_new)*100)
end do
close (10,status='keep')

write(9,*) 'Last ocean index:', ip1
write(9,*) 'Last land index:', ip2
write(9,*) 'NTILES LAND:', nland_new
write(*,*)

allocate ( lats_new(nland_new) )
allocate ( lons_new(nland_new) )
allocate (   fr_new(nland_new) )
allocate (   ii_new(nland_new) )
allocate (   jj_new(nland_new) )

lats_new = pack(lats_tmp, mask=typ_tmp .eq. 100)
lons_new = pack(lons_tmp, mask=typ_tmp .eq. 100)
  fr_new = pack(  fr_tmp, mask=typ_tmp .eq. 100)
  ii_new = pack(  ii_tmp, mask=typ_tmp .eq. 100)
  jj_new = pack(  jj_tmp, mask=typ_tmp .eq. 100)

write(9,*) 'lats', size(lats_new), minval(lats_new), maxval(lats_new)
write(9,*) 'lons', size(lons_new), minval(lons_new), maxval(lons_new)
write(9,*) 'fr  ',   size(fr_new), minval(fr_new), maxval(fr_new) 
write(9,*) 'ii  ',   size(ii_new), minval(ii_new), maxval(ii_new)
write(9,*) 'jj  ',   size(jj_new), minval(jj_new), maxval(jj_new)

deallocate(lats_tmp)
deallocate(lons_tmp)
deallocate(  fr_tmp)
deallocate(  ii_tmp)
deallocate(  jj_tmp)
deallocate( typ_tmp)

! -------------------------------------------------------------------------------
! 3. Read in the old restart from a previous run
!    Here, I separate the parameters and prognostic variables.  Some of the 
!    prognostic variables are printed out by catch-finalize as var(ntiles, 4)
!    This routine takes that into account, and I put the parameters and 
!    prognostics in separate arrays.  Then, the parameters will be replaced by
!    something Sarith makes, while the prognostics will be regridded.  If you 
!    wish, you can also retain the soil parameters in the trivial case (eg.
!    you want to keep same land specification but just adjust the initialization
!    to a different date) 
! -------------------------------------------------------------------------------

allocate(oldparmvars(nland_old,        numparmrecs))
allocate(oldprogvars(nland_old,allrecs-numparmrecs))

allocate( oldallvars(nland_old,allrecs))
allocate(tiletilevar(nland_old, numsubtiles))
allocate(    tilevar(nland_old))

open(unit=30, file=trim(old_restart),form='unformatted')

write (*,*)
write (*,*) '---------------------------------------------------------------------'
write (*,*) 'Reading old restart from:'
write (*,*) trim(old_restart)

write (9,*)
write (9,*) '---------------------------------------------------------------------'
write (9,*) 'Reading '//trim(old_restart)
write (9,*) 'Sizes', size(tiletilevar), size(tilevar)

 write(*, 70, advance=trim(flag))  
open(unit=65, file='old_catch.dat' ,form='unformatted')

nta=1
do n=1, numrecs   
   if (ttlookup(n)) then 
      read(30) tiletilevar
      do nn=1, numsubtiles
         oldallvars(:,nta)=tiletilevar(:,nn)
         write (65)        tiletilevar(:,nn)  ! Write Grads-Formatted Catchment File
         nta=nta+1
      enddo 
   else
      read(30) tilevar
      oldallvars(:,nta)=tilevar(:)
      write (65)        tilevar(:)            ! Write Grads-Formatted Catchment File
      nta=nta+1
   endif
   write(*, 50, advance=trim(flag)) bak, floor(float(n)/float(numrecs)*100)
enddo  

close(30)
deallocate(tiletilevar)
deallocate(tilevar)

write(*,*) 'Separating parameter and prognostic variables'
write(9,*) 'Separating parameter and prognostic variables'

do n=1, numparmrecs
   oldparmvars(:,n)=oldallvars(:,n)
enddo
do n=1, allrecs-numparmrecs
   oldprogvars(:,n)=oldallvars(:,n+numparmrecs)
end do 

loc = 0
do n=1,numrecs
                      nta = 1
    if( ttlookup(n) ) nta = numsubtiles
    do nn = 1,nta
      loc = loc+1
      if( loc.le.numparmrecs ) then
          write(9,*) '  Transferred old parameter  (',n,',',nn,') ',  &
	                minval(oldallvars(:,loc)), maxval(oldallvars(:,loc)) 
      else
          write(9,*) '  Transferred old prognostic (',n,',',nn,') ',  &
	                minval(oldallvars(:,loc)), maxval(oldallvars(:,loc)) 
      endif
    enddo
enddo


! -------------------------------------------------------------------------------
! 4. Read in the soil parameter variables (there are 29 of them) from Sarith
! vegetation type is also read in here, from an old Aries format (this needs
! to be changed, so vegtype is in .til file!)
! -------------------------------------------------------------------------------

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

write(*,*)
write(*,*) '---------------------------------------------------------------------'
write(9,*) 
write(9,*) '---------------------------------------------------------------------'
write(9,*) 'Reading Sarith parameters from:'
write(9,*) trim(newtildir)
write(9,*) 'Sample output ... '
write(*,*) 'Reading Sarith parameters from:'
write(*,*) trim(newtildir)
write(*,*) 'Sample output ... '

open(unit=21, file=trim(newtildir) // '/' //'mosaic_veg_typs_fracs',form='formatted')
open(unit=22, file=trim(newtildir) // '/' //'bf.dat'               ,form='formatted')
open(unit=23, file=trim(newtildir) // '/' //'soil_param.dat'       ,form='formatted')
open(unit=24, file=trim(newtildir) // '/' //'ar.new'               ,form='formatted')
open(unit=25, file=trim(newtildir) // '/' //'ts.dat'               ,form='formatted')
open(unit=26, file=trim(newtildir) // '/' //'tau_param.dat'        ,form='formatted')

 write(*, 80, advance=trim(flag))  

do n=1,nland_new 
!  read (21, *) catindex21, catid, ity_int(n), ity2, frc1, frc2, rdum
   read (21, *) catindex21, catid, ity_int(n), ity2, frc1, frc2  ! version 2 doesnt have rdum variable
   ITY0(n)=1.0*ity_int(n)
   read (22, *) catindex22, catid, GNU(n), BF1(n), BF2(n), BF3(n)
   read (23, *) catindex23, catid, idum, idum, BEE(n), PSIS(n), POROS(n), COND(n), WPWET(n), DP2BR(n)

   read (24, *) catindex24, catid, rdum, ARS1(n), ARS2(n), ARS3(n),          &
                                         ARA1(n), ARA2(n), ARA3(n), ARA4(n), &
					 ARW1(n), ARW2(n), ARW3(n), ARW4(n)

   read (25, *) catindex25, catid, rdum, TSA1(n), TSA2(n), TSB1(n), TSB2(n)

   read (26, *) catindex26, catid, ATAU2(n), BTAU2(n), rdum, rdum

   checksum=catindex21+catindex22+catindex23+catindex24+catindex25+catindex26-6*(n+ip1)
   zdep2=1000.
   zdep3=amax1(1000.,DP2BR(n))
   if (zdep2 .gt.0.75*zdep3) then
       zdep2  =  0.75*zdep3              
   end if
   zdep1=20.
   zmet=zdep3/1000.
   term1=-1.+((PSIS(n)-zmet)/PSIS(n))**((BEE(n)-1.)/BEE(n))
   term2=PSIS(n)*BEE(n)/(BEE(n)-1)
   VGWMAX(n)=POROS(n)*zdep2   
   CDCR1(n)=1000.*POROS(n)*(zmet-(-term2*term1))   
   CDCR2(n)=(1.-WPWET(n))*POROS(n)*zdep3
   if (checksum .ne. 0) then 
      write(9,*) 'Catchment id mismatch with following id list at n=', n
      write(9,*) catindex22, catindex23, catindex24, catindex25, catindex26, ip1+n
      write(*,*) 'Halted on catchment mismatch'
      STOP
   else 
      if (modulo(n, 1000).eq.1 .or. n.eq.qtile ) then
         write(9,*)
         write(9,*) n, 'mosaic_vegtype: ', ity_int(n)
         write(9,*) n, 'bf.dat:         ', catindex22, catid, GNU(n), BF1(n), BF2(n)
         write(9,*) n, 'bf.dat:         ', BF3(n)
         write(9,*) n, 'soil_param.dat: ', catindex23, catid, rdum, BEE(n), PSIS(n)
         write(9,*) n, 'soil_param.dat: ', POROS(n), COND(n), WPWET(n), DP2BR(n)
         write(9,*) n, 'ar.dat:         ',  catindex24, catid, rdum, ARS1(n), ARS2(n)
         write(9,*) n, 'ar.dat:         ', ARS3(n), ARA1(n), ARA2(n), ARA3(n), ARA4(n)
         write(9,*) n, 'ar.dat:         ', ARW1(n), ARW2(n), ARW3(n), ARW4(n)
         write(9,*) n, 'ts.dat:         ', catindex25, catid, rdum, TSA1(n), TSA2(n)
         write(9,*) n, 'ts.dat:         ', TSB1(n), TSB2(n)
         write(9,*) n, 'tau_param.dat:  ', catindex26, catid, ATAU2(n), BTAU2(n)
         write(9,*) n, 'Computed:       ', VGWMAX(n), CDCR1(n), CDCR2(n) 
      end if
   endif
   write(*, 50, advance=trim(flag)) bak, floor(float(n)/float(nland_new)*100)
end do

close (21)
close (22)
close (23)
close (24)
close (25)
close (26)

! -------------------------------------------------------------------------------
! 5. Regrid all  variables to im_gcm_oldXjm_gcm_old grid
!    Then, find interpolated values based upon the centroids of tiles in
!    new .til definitions.  Missing values:  if a single tile is missing,
!    it's influence on the gridded value is ignored, except if there are no
!    non-missing values in an i,j cell, then the new tile is defined as missing
!   
!    Alternatives for future development:
!
!    a. Nearest neighbor
!    b. Nearest neighbor of same/similar vegetation type, latitude, etc.
!    c. Gridding, ungridding (this is done currently)
!    d. Krieging of some kind, pick a radius of influence and weigh by inverse
!       square of distance, or limit to veg type, or whatever.  
!  
! -------------------------------------------------------------------------------

open(unit=8, file=trim(old_diag_grids),form='unformatted')

allocate(    tmp_sum(im_gcm_old,jm_gcm_old))
allocate(    tmp_wgt(im_gcm_old,jm_gcm_old))
allocate(oldvargrids(im_gcm_old,jm_gcm_old,allrecs)) 
allocate(    dumgrid(im_gcm_old,jm_gcm_old)) 

loc = 0
do v=1, numrecs
                      nta = 1
    if( ttlookup(v) ) nta = numsubtiles
    do nn = 1,nta
      loc = loc+1
      tmp_sum(:,:)=0.0
      tmp_wgt(:,:)=0.0
      do n=1, nland_old
         val0=oldallvars(n,loc)
         ii0=ii_old(n)
         jj0=jj_old(n)
         fr0=fr_old(n)
         if (abs(val0-ESMF_MISSING) .gt. tol) then 
            tmp_sum(ii0,jj0) = tmp_sum(ii0,jj0) + fr0*val0
            tmp_wgt(ii0,jj0) = tmp_wgt(ii0,jj0) + fr0
         else
            print *, 'Old_Catch_Val = ',val0,' n = ',n,' loc = ',loc
         endif      
      enddo
      do j=1,jm_gcm_old
      do i=1,im_gcm_old
         if (tmp_wgt(i,j) .gt. tol) then 
             oldvargrids(i,j,loc)=tmp_sum(i,j)/tmp_wgt(i,j)
         else 
             oldvargrids(i,j,loc)=ESMF_MISSING
         endif
                 dumgrid(i,j)    =oldvargrids(i,j,loc)
      enddo
      enddo
    write (8)    dumgrid
    enddo
enddo 
   
deallocate(dumgrid)
deallocate(tmp_sum)
deallocate(tmp_wgt)
close(8)

allocate( corners_lookup(nland_new, 4) )
allocate( weights_lookup(nland_new, 4) )

lat00_old = -90
   dx_old = (360.0)/ im_gcm_old
   dy_old = (180.0)/(jm_gcm_old-1)

if (old_dateline .eq. 'DC') then
   lon00_old = -180
   lonIM_old =  180-dx_old
else 
   lon00_old = -180+0.5*dx_old
   lonIM_old =  180-0.5*dx_old
end if

write(*,*)
write(*,*) '---------------------------------------------------------------------'
write(*,*) 'Computing interpolation lookup table for new tiles'
write(9,*)
write(9,*) '---------------------------------------------------------------------'
write(9,*) 'Computing interpolation lookup table for new tiles'
 write(*, 90, advance=trim(flag))  

do n=1,nland_new
   lat0=lats_new(n)                                !  latitude of tile centroid to find
   lon0=lons_new(n)                                !  longitude of tile centroid
   if ((lon0 .gt. lonIM_old) .or. (lon0 .lt. lon00_old)) then 
      ia=im_gcm_old                                           
      ib=1                               
      lona=lonIM_old                
      lonb=lon00_old
   else 
      ia=floor((lon0-lon00_old)/dx_old)+1
      ib=ia+1
      lona=(ia-1)*dx_old+lon00_old  
      lonb=lona+dx_old
   end if
   ja=floor((lat0-lat00_old)/dy_old)+1          !  left bottom corner y coordinate
   jb=ja+1                                      !  right top corner y coordinate
   lata=(ja-1)*dy_old+lat00_old                 !  latitude of left bottom corner
   latb=lata+dy_old   

   if( ia.lt.1 .or. ia.gt.im_gcm_old .or. &
       ib.lt.1 .or. ib.gt.im_gcm_old .or. &
       ja.lt.1 .or. ja.gt.jm_gcm_old .or. &
       jb.lt.1 .or. jb.gt.jm_gcm_old ) then
       print *, 'Warning, bad indicies!'
       print *, 'New Land variable: ',n,ia,ib,ja,jb
       stop
   endif

   if (modulo(n, 1000).eq.1 .or. n.eq.qtile) then
      write (9,*)
      write (9,*) n, lona, lon0, lonb, ia, ib
      write (9,*) n, lata, lat0, latb, ja, jb
   end if
   corners_lookup(n,1)=ia
   corners_lookup(n,2)=ib
   corners_lookup(n,3)=ja
   corners_lookup(n,4)=jb
   waa=sqrt((lat0-lata)**2+(lon0-lona)**2)
   wab=sqrt((lat0-lata)**2+(lon0-lonb)**2)
   wba=sqrt((lat0-latb)**2+(lon0-lona)**2)
   wbb=sqrt((lat0-latb)**2+(lon0-lonb)**2)
   wsum=waa+wab+wba+wbb
   weights_lookup(n,1)=waa
   weights_lookup(n,2)=wab
   weights_lookup(n,3)=wba
   weights_lookup(n,4)=wbb
   write(*, 50, advance=trim(flag)) bak, floor(float(n)/float(nland_new)*100)
end do 

allocate(newallvars(nland_new, allrecs)) 
! new allvars allocated by number of new land tiles X number of total restart records
write(*,*)
write(*,*) '---------------------------------------------------------------------'
write(*,*) 'Interpolating prognostic records to new tile definitions'
write(9,*)
write(9,*) '---------------------------------------------------------------------'
write(9,*) 'Interpolating prognostic records to new tile definitions'
 write(*, 90, advance=trim(flag)) 
 write(*, 50, advance=trim(flag)) bak, 0

do v=31, allrecs
   do n=1, nland_new
       ia=corners_lookup(n,1)
       ib=corners_lookup(n,2)
       ja=corners_lookup(n,3)
       jb=corners_lookup(n,4)
      waa=weights_lookup(n,1)
      wbb=weights_lookup(n,4)
      wab=weights_lookup(n,2)
      wba=weights_lookup(n,3)
      wsum=0
      tempval=0
      vaa=oldvargrids(ia,ja, v)
      vbb=oldvargrids(ib,jb, v)
      vab=oldvargrids(ia,jb, v)
      vba=oldvargrids(ib,ja, v)
      if (abs(vaa-ESMF_MISSING) .gt. tol) then
         tempval=tempval+vaa*waa
         wsum=waa+wsum
      end if
      if (abs(vab-ESMF_MISSING) .gt. tol) then
         tempval=tempval+vab*wab
         wsum=wab+wsum
      end if
      if (abs(vba-ESMF_MISSING) .gt. tol) then
         tempval=tempval+vba*wba
         wsum=wsum+wba
      end if
      if (abs(vbb-ESMF_MISSING) .gt. tol) then
         tempval=tempval+vbb*wbb
         wsum=wsum+wbb
      end if
      if (abs(wsum) .lt. tol) then 
              dmin =  1e15
              nmin =  0
              do nn = 1,nland_old
              dist = sqrt(  (lats_old(nn)-lats_new(n))**2 &
                          + (lons_old(nn)-lons_new(n))**2 )
              if( dist.lt.dmin ) then
                  nmin = nn
                  dmin = dist
              endif
              enddo
         tempval=oldallvars(nmin,v)  ! Find nearest old tile to new tile
         print *, 'NewVal = ',tempval,' nmin = ',nmin,' loc = ',v
         print *, 'newlat = ',lats_new(n),' oldlat = ',lats_old(nmin)
         print *, 'newlon = ',lons_new(n),' oldlon = ',lons_old(nmin)
         print *
      else 
         tempval=tempval/(1.0*wsum)
      end if
      newallvars(n, v)=tempval 
      if (modulo(n, 1000) .eq. 1) then
         write(9,*)
         write(9,*) n, 'Interpolation summary'
         write(9,*) n, 'Weights:', waa, wab, wba, wbb
         write(9,*) n, 'Values:',  vaa, vab, vba, vbb
         write(9,*) n, 'Results:', tempval, wsum
      end if
   end do

   write(*,*) v
   write(*, 50, advance=trim(flag)) bak, floor((float(v-31)/float(allrecs-31)*100))
end do 

! I am now finished with the old tiles, get rid of them
deallocate(oldprogvars, oldparmvars, oldallvars)
deallocate(corners_lookup, weights_lookup)

! -------------------------------------------------------------------------------
! 6. Create the restart from stored values
!    a. 29 Sarith tilespace records from his parameter files
!    b. The vegetation type from Sarith's mosaic_veg_typ_file 
!       (I have used the PRIMARY veg type for this work, as opposed 
!       to the second one that also has a fraction. I am assuming that 
!       the catchment fraction is totally composed of the PRIMARY veg type)
!    c. The modified/regridded prognostic variables that have been regridded
!       At this point, the variable newallvars contains ALL records for the
!       restart, including estimates of Sarith's parameters based upon 
!       the old values interpolated from the old restart.  These are skipped, but
!       might be useful for comparison in a debugging situation
! -------------------------------------------------------------------------------

open(unit=41, file=trim(new_restart),form='unformatted')
open(unit=66, file='new_catch.dat'  ,form='unformatted')

! replace the old interpolated parameters in newallvars with the new Sarith ones

   write(9,*) '  Min/Max for ARS1: ', minval(ARS1), maxval(ARS1) 
   write(9,*) '  Min/Max for ARS2: ', minval(ARS2), maxval(ARS2) 
   write(9,*) '  Min/Max for ARS3: ', minval(ARS3), maxval(ARS3) 

newallvars(:,1)=BF1
newallvars(:,2)=BF2
newallvars(:,3)=BF3
newallvars(:,4)=VGWMAX
newallvars(:,5)=CDCR1
newallvars(:,6)=CDCR2
newallvars(:,7)=PSIS
newallvars(:,8)=BEE
newallvars(:,9)=POROS 
newallvars(:,10)=WPWET
newallvars(:,11)=COND
newallvars(:,12)=GNU
newallvars(:,13)=ARS1
newallvars(:,14)=ARS2
newallvars(:,15)=ARS3
newallvars(:,16)=ARA1
newallvars(:,17)=ARA2
newallvars(:,18)=ARA3
newallvars(:,19)=ARA4
newallvars(:,20)=ARW1
newallvars(:,21)=ARW2
newallvars(:,22)=ARW3
newallvars(:,23)=ARW4
newallvars(:,24)=TSA1
newallvars(:,25)=TSA2
newallvars(:,26)=TSB1
newallvars(:,27)=TSB2
newallvars(:,28)=ATAU2
newallvars(:,29)=BTAU2
newallvars(:,30)=ITY0

write(*,*)
write(*,*) '---------------------------------------------------------------------'
write(*,*) 'Writing new restart'
write(9,*)
write(9,*) '---------------------------------------------------------------------'
write(9,*) 'Writing new restart'

loc = 0
do v=1,numrecs
                      nta = 1
    if( ttlookup(v) ) nta = numsubtiles
    allocate( dumtile(nland_new,nta) )
    do nn = 1,nta
      loc = loc+1
      dumtile(:,nn) = newallvars(:,loc)
      write (66)         dumtile(:,nn)  ! Write Grads-Formatted Catchment File
    enddo
    write (41) dumtile
    write(9,*) 'NEW RESTART RECORD #', v , ' Size = ',size(dumtile)
    deallocate ( dumtile )
enddo 
close(41)

loc = 0
do n=1,numrecs
                      nta = 1
    if( ttlookup(n) ) nta = numsubtiles
    do nn = 1,nta
      loc = loc+1
      if( loc.le.numparmrecs ) then
          write(9,*) '  Transferred new parameter  (',n,',',nn,') ',  &
	                minval(newallvars(:,loc)), maxval(newallvars(:,loc)) 
      else
          write(9,*) '  Transferred new prognostic (',n,',',nn,') ',  &
	                minval(newallvars(:,loc)), maxval(newallvars(:,loc)) 
      endif
    enddo
enddo


! -------------------------------------------------------------------------------
! 7. Save a gridded copy of the new restart on rectangular grid found in the 
!    new .til file definitions.  This can be used to check the results.
! -------------------------------------------------------------------------------

open(unit=42, file=trim(new_diag_grids),form='unformatted')

allocate(   tmp_sum(im_gcm_new, jm_gcm_new))
allocate(   tmp_wgt(im_gcm_new, jm_gcm_new))
allocate(   dumgrid(im_gcm_new, jm_gcm_new)) 

loc = 0
do v=1, numrecs
                      nta = 1
    if( ttlookup(v) ) nta = numsubtiles
    do nn = 1,nta
      loc = loc+1
     tmp_sum(:,:)=0.0
     tmp_wgt(:,:)=0.0
     do n=1, nland_new
        val0=newallvars(n,loc)
        ii0=ii_new(n)
        jj0=jj_new(n)
        fr0=fr_new(n)
        if (abs(val0-ESMF_MISSING) .gt. tol) then 
           tmp_sum(ii0,jj0)=tmp_sum(ii0,jj0)+fr0*val0
           tmp_wgt(ii0,jj0)=tmp_wgt(ii0,jj0)+fr0
        endif      
     enddo
     do i=1,im_gcm_new
     do j=1,jm_gcm_new
        if (tmp_wgt(i,j) .gt. tol) then 
            dumgrid(i,j)=tmp_sum(i,j)/tmp_wgt(i,j)
        else 
            dumgrid(i,j)=ESMF_MISSING
        endif
     enddo
     enddo
    write (42)  dumgrid
    enddo
enddo 
close(42)

deallocate(tmp_sum)
deallocate(tmp_wgt)
deallocate( dumgrid ) 

deallocate(BF1,   BF2,   BF3,  VGWMAX)
deallocate(CDCR1, CDCR2, PSIS, BEE) 
deallocate(POROS, WPWET, COND, GNU)
deallocate(ARS1,  ARS2,  ARS3)
deallocate(ARA1,  ARA2,  ARA3)
deallocate(ARA4,  ARW1,  ARW2, ARW3, ARW4)
deallocate(TSA1,  TSA2,  TSB1, TSB2)
deallocate(DP2BR, ATAU2, BTAU2)
deallocate(ITY0, ity_int)

deallocate(lats_old)
deallocate(lons_old)
deallocate(fr_old)
deallocate(ii_old)
deallocate(jj_old)
deallocate(lats_new)
deallocate(lons_new)
deallocate(fr_new)
deallocate(ii_new)
deallocate(jj_new)
deallocate(ttlookup)

40 FORMAT(' Percent tile definitions read:      ')
50 FORMAT(A4, I3.3, '%') 
60 FORMAT(' Percent MODIS data read:      ')
70 FORMAT(' Percent restart read:      ')
80 FORMAT(' Percent Sarith catchment parameters read:      ')
90 FORMAT(' Percent completed:      ')
END
