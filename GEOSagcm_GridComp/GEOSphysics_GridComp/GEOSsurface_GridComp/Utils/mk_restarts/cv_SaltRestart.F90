

program cv_SaltRestart

  implicit none


  integer, parameter :: char_len  = 80, &
                        char_len_long  = 256, &
                        log_kind  = kind(.true.), &
                        int_kind  = selected_int_kind(6), &
                        real_kind = selected_real_kind(6), &
                        dbl_kind  = selected_real_kind(13), &
                        r16_kind  = selected_real_kind(26)


  integer :: i, n, k, nxt, argc, ntiles, nrec
  character*128 :: Outfile, InFile1, InFile2, length
  character*128 :: tilefile
  character*1            :: Opt
  character*128          :: arg


  real*4, allocatable, dimension(:  ) :: HW,TW,SW,HI,TI,SI
  real*4, allocatable, dimension(:  ) :: rainp, rainn, snowp, snown, rrp, rrn
  real*4, allocatable, dimension(:  ) :: swradp, swradn, lwradp, lwradn, t10p, t10n
  real*4, allocatable, dimension(:,:) :: tauage
  real*4, allocatable, dimension(:)   :: slmask 
  real*4, allocatable, dimension(:,:) :: QS,CH,CM,CQ,Z0,WW
  real*4, allocatable, dimension(:  ) :: TWMTS, DTWARM 
  real*4, allocatable, dimension(:)   :: sst, lons, lats 
  real*4, allocatable, dimension(:,:) :: Tsc, Ts
  real*4, allocatable, dimension(:)   :: frice
  real*4, allocatable, dimension(:    ) :: X1D4

  real*8, allocatable, dimension(:    ) :: X1D8
  real*8, allocatable, dimension(:,  :) :: & 
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon , & ! volume per unit area of snow         (m)
         trcrn , & ! ice tracers
                   ! 1: surface temperature of ice/snow (C)
         volpondn, &
         apondn, & 
         hpondn, & 
         eicen , & ! energy of melting for each ice layer  (J/m^2)
         esnon     ! energy of melting for each ice layer  (J/m^2)
  real*4                                 :: lono, lato

  character*128 :: Usage1="Usage: cv_SaltRestart -A -s MERRA2-SaltInternal MERRA2-SaltImport -t tilefile"
  character*128 :: Usage2="Usage: cv_SaltRestart -C -s Coupled-SaltInternal -t tilefile"
  integer, parameter   :: ncat = 5
  integer, parameter   :: nsub = ncat+1
  integer, parameter   :: nilyr = 4
  integer, parameter   :: nslyr = 1
  integer, parameter   :: ntrcr = 1
  integer, parameter   :: ntilyr = nilyr*ncat
  integer, parameter   :: ntslyr = nslyr*ncat
  real,    parameter   :: Tffresh   = 273.15 ! freezing temp of fresh ice (K)
  integer ::  ilyr1 (ncat), & ! array position of top ice layer in each cat
              ilyrn (ncat), & ! array position of bottom ice layer in each cat
              slyr1 (ncat), & ! array position of top snow layer in each cat
              slyrn (ncat)    ! array position of bottom snow layer in each cat


  logical               :: source_is_coupled


  source_is_coupled = .false.

  I = iargc()

  print*, '# of arguments : ', I

  if(I < 5 .or. I > 6) then
     print *, trim(Usage1)
     print *, ' ========= or ==========='
     print *, trim(Usage2)
     call exit(1)
  endif

  nxt = 1

  do while(nxt<=I)
       call getarg(nxt,arg)
       if(arg(1:1)/='-') then
          print*, 'first letter /= -'
          print *, trim(Usage1)
          print *, ' ========= or ==========='
          print *, trim(Usage2)
          call exit(1)
       endif
       opt=arg(2:2)
       if(len(trim(arg))==2) then
          if(scan(opt,'AC')==0) then
             nxt = nxt + 1
             call getarg(nxt,arg)
          end if
       else
          arg = arg(3:)
       end if
       select case (opt)
       case ('A')
          source_is_coupled = .false.
       case ('C')
          source_is_coupled = .true.
       case ('s')
          InFile1 = arg
          if(.not. source_is_coupled) then
             nxt = nxt + 1
             call getarg(nxt,arg)
             InFile2 = arg
          endif  
       case ('t')
          TileFile = arg
       case default
          print*, 'unrecognized letter'
          print *, trim(Usage1)
          print *, ' ========= or ==========='
          print *, trim(Usage2)
          call exit(1)
       end select

       nxt = nxt + 1
  enddo

  if(source_is_coupled) then
      print*, 'source is from coupled model'
  else
      print*, 'source is from AMIP'
  endif
  print*, trim(InFile1)
  if(.not. source_is_coupled) print*, trim(InFile2)
  print*, trim(TileFile)






!#if 0


  if(.not. source_is_coupled)  then
    print*, 'conversion from AMIP style restarts .NOT. implemented yet!!'
    call exit(1)
  else

  ilyr1(1) = 1                       ! if nilyr  = 4
  ilyrn(1) = nilyr
  do n = 2, ncat                     !   ilyrn = { 4,8,12} etc
     ilyr1(n) = ilyrn(n-1) + 1
     ilyrn(n) = ilyrn(n-1) + nilyr
  enddo
  slyr1(1) = 1
  slyrn(1) = nslyr
  do n = 2, ncat
     slyr1(n) = slyrn(n-1) + 1
     slyrn(n) = slyrn(n-1) + nslyr
  enddo


  print*,(ilyr1(n), n=1,ncat)


  ntiles = GetNumTiles(tilefile)


  print *, "Processing restarts with ", ntiles, "tiles" 

  allocate(HW(ntiles),TW(ntiles),SW(ntiles),HI(ntiles),TI(ntiles),SI(ntiles), &
           QS(ntiles,nsub),CH(ntiles,nsub),CM(ntiles,nsub), &
           CQ(ntiles,nsub),Z0(ntiles,nsub),WW(ntiles,nsub)  )

  allocate(sst(ntiles), lons(ntiles), lats(ntiles),  frice(ntiles))
  allocate(TWMTS(ntiles), DTWARM(ntiles))
  

  allocate(X1D8(ntiles))
  allocate(X1D4(ntiles))
  allocate(aicen(ntiles,nsub), vicen(ntiles,ncat), vsnon(ntiles,ncat), trcrn(ntiles,ncat), &
           eicen(ntiles,ncat*nilyr), esnon(ntiles,ncat*nslyr))
  allocate(volpondn(ntiles,ncat), apondn(ntiles,ncat), hpondn(ntiles,ncat))
  allocate(Tsc(ntiles,nsub))
  allocate(Ts(ntiles,nsub))
  allocate(tauage(ntiles,ncat), slmask(ntiles))

  call GetTileLonLats(tilefile, ntiles, lons, lats) 

  print*, "InFile1: ",InFile1 
  print*, "TileFile: ",tilefile

  i = index(InFile1,'/',back=.true.)

  open(10, file=InFile1, form="unformatted",status='old', convert='little_endian')
  !open(30, file=InFile2, form="unformatted",status='old', convert='little_endian')
  open(20, file="OutData/"//trim(InFile1(i+1:)), form="unformatted", &
          status='unknown', convert='little_endian')




  read(10)  hi
  read(10)  ti
  read(10)  si
  read(10)  hw
  read(10)  sw
  read(10)  tw

  nrec = 6

  ! tskin
  read(10)  X1D4
  Ts(:,1) = X1D4
  do i=1,ncat
     read(10)  X1D4
     Ts(:,i+1) = X1D4
  enddo

! tskinc
  read(10)  X1D4
  Tsc(:,1) = X1D4
  do i=1,ncat
     read(10)  X1D4
     Tsc(:,i+1) = X1D4
  enddo

! FR

  do i=1,nsub
     read(10)  X1D8
     aicen(:,i) = X1D8
  enddo

  do i=1,ncat
     read(10)  X1D8
     vicen(:,i) = X1D8
  enddo

  do i=1,ncat
     read(10)  X1D8
     vsnon(:,i) = X1D8
  enddo

  do i=1,ncat
     read(10)  X1D8
     volpondn(:,i) = X1D8
  enddo

  do i=1,ncat
     read(10)  X1D8
     apondn(:,i) = X1D8
  enddo

  do i=1,ncat
     read(10)  X1D8
     hpondn(:,i) = X1D8
  enddo

  do k=1,nilyr
       do i=1,ncat
         read(10)  X1D8
         eicen(:,ilyr1(i)+k-1) = X1D8
       enddo
  enddo

  do k=1,nslyr
     do i=1,ncat
       read(10)  X1D8
       esnon(:,slyr1(i)+k-1) = X1D8
     enddo
  enddo


!ERGSUM
  do i=1,ncat
     read(10)  X1D4
     nrec = nrec + 1
  enddo

!TAUAGE
  do i=1,ncat
     read(10)  X1D4
     tauage(:,i) = X1D4
     nrec = nrec + 1
  enddo


!FRZMLT
  read(10)  X1D4
  nrec = nrec + 1

!QS
  read(10)  X1D4
  QS(:,1)  = X1D4
  nrec = nrec + 1
  do i=1,ncat
     read(10)  X1D4
     QS(:,i+1)  = X1D4
     nrec = nrec + 1
  enddo

!CH
  read(10)  X1D4
  CH(:,1)  = X1D4
  nrec = nrec + 1
  do i=1,ncat
     read(10)  X1D4
     CH(:,i+1)  = X1D4
     nrec = nrec + 1
  enddo


!CM
  read(10)  X1D4
  CM(:,1)  = X1D4
  nrec = nrec + 1
  do i=1,ncat
     read(10)  X1D4
     CM(:,i+1)  = X1D4
     nrec = nrec + 1
  enddo

!CQ
  read(10)  X1D4
  CQ(:,1)  = X1D4
  nrec = nrec + 1
  do i=1,ncat
     read(10)  X1D4
     CQ(:,i+1)  = X1D4
     nrec = nrec + 1
  enddo

!Z0
  read(10)  X1D4
  Z0(:,1)  = X1D4
  nrec = nrec + 1
  do i=1,ncat
     read(10)  X1D4
     Z0(:,i+1)  = X1D4
     nrec = nrec + 1
  enddo

!WW
  read(10)  X1D4
  WW(:,1)  = X1D4
  nrec = nrec + 1
  do i=1,ncat
     read(10)  X1D4
     WW(:,i+1)  = X1D4
     nrec = nrec + 1
  enddo

!HFLUX
  read(10)  X1D4
  nrec = nrec + 1
!WATERFLUX
  read(10)  X1D4
  nrec = nrec + 1
!SALTFLUX
  read(10)  X1D4
  nrec = nrec + 1
!SLMASK
  read(10)  X1D4
  slmask = X1D4
  nrec = nrec + 1

  print*, 'Finished reading old saltwater internal restart file'

  print*, 'Start writing new CICEThermo internal restart file ...' 

  write(20) hw
  write(20) tw
  write(20) sw
  write(20) hi
  nrec = 4
  do i=1,ncat
     X1D4 = Ts(:,i+1)
     write(20) X1D4
     nrec = nrec + 1
  enddo
  write(20) si
  nrec = nrec + 1

!QS
  !do i=1,nsub
  !   write(20) QS(:,i)
  !   nrec = nrec + 1
  !enddo
  ! for tiletile vars, they need to be written out as one single record
  write(20) QS

!CH
  !do i=1,nsub
  !   write(20) CH(:,i)
  !   nrec = nrec + 1
  !enddo
  write(20) CH

!CM
  !do i=1,nsub
  !   write(20) CM(:,i)
  !   nrec = nrec + 1
  !enddo
  write(20) CM

!CQ
  !do i=1,nsub
  !   write(20) CQ(:,i)
  !   nrec = nrec + 1
  !enddo
  write(20) CQ

!Z0
  !do i=1,nsub
  !   write(20) Z0(:,i)
  !   nrec = nrec + 1
  !enddo
  write(20) Z0

!WW
  !do i=1,nsub
  !   write(20) WW(:,i)
  !   nrec = nrec + 1
  !enddo
  write(20) WW

!TWMTS
  X1D4 = 0.0
  write(20) X1D4
  nrec = nrec + 1
 

! FR
  do i=1,nsub
     write(20) aicen(:,i) 
     nrec = nrec + 1
  enddo

! VOLICE
  do i=1,ncat
     write(20) vicen(:,i) 
     nrec = nrec + 1
  enddo

! VOLSNO
  do i=1,ncat
     write(20) vsnon(:,i) 
     nrec = nrec + 1
  enddo

! VOLPOND
  do i=1,ncat
     write(20) volpondn(:,i) 
     nrec = nrec + 1
  enddo

! APOND
  do i=1,ncat
     write(20) apondn(:,i) 
     nrec = nrec + 1
  enddo

! HPOND
  do i=1,ncat
     write(20) hpondn(:,i) 
     nrec = nrec + 1
  enddo

!ERGICE
  do k=1,nilyr
       do i=1,ncat
         X1D8 = eicen(:,ilyr1(i)+k-1) 
         write(20) X1D8 
         nrec = nrec + 1
       enddo 
  enddo

!ERGSNO
  do k=1,nslyr
     do i=1,ncat
       X1D8 = esnon(:,slyr1(i)+k-1)
       write(20) X1D8 
       nrec = nrec + 1
     enddo
  enddo

!TAUAGE
  do i=1,ncat
     write(20) tauage(:,i) 
     nrec = nrec + 1
  enddo

!SLMASK
  write(20) slmask 
  nrec = nrec + 1

  close(10)
  CLOSE(20)



  print *, 'Wrote ', nrec, ' records'

 lono = -164.2552
 lato = 59.18839

  do k=1,ntiles
      if(abs(lons(k)-lono)<1.e-4 .and. abs(lats(k)-lato)< 1.e-4) then
         print*,(aicen(k,n), n=1,ncat)
         print*,(vicen(k,n), n=1,ncat)
         print*,frice(k),   tw(k),   ti(k)
         do i=1,ncat
            print*,(eicen(k,ilyr1(i)+n-1),n=1,nilyr)
         enddo 
      endif
  enddo
  !deallocate(rainp, rainn, snowp, snown, rrp, rrn,       &
  !         swradp, swradn, lwradp, lwradn, t10p, t10n, &
  !         q10p, q10n,  u10p, u10n,  v10p, v10n,       &
  !         slpp, slpn,  sssp, sssn )
  deallocate(sst, lons, lats, frice)
  deallocate(aicen, vicen, vsnon, trcrn, &
           eicen, esnon)
  deallocate(Tsc)

  endif 

!#endif

contains



integer function GetNumTiles(tilefile)
  implicit none
  character*128, intent(in) :: tilefile
! local
  integer                  :: ntiles
  integer                  :: mark
  integer                  :: dum, n, nt
  character*128            :: dumstr
  
  
  open(10, file=tilefile, form="formatted",status='old')

  read(10, fmt=*) nt
  read(10, fmt=*) dum
  read(10, fmt=*) dumstr
  read(10, fmt=*) dum
  read(10, fmt=*) dum
  read(10, fmt=*) dumstr
  read(10, fmt=*) dum
  read(10, fmt=*) dum
  ntiles = 0
  do n=1,nt
    read(10, fmt=*) mark
    if(mark .eq. 0) ntiles = ntiles + 1
  enddo
  close(10)
 
  GetNumTiles = ntiles 

end function GetNumTiles

subroutine GetTileLonLats(tilefile, ntiles, Lons, Lats)
  implicit none
  character*128, intent(in) :: tilefile
  integer, intent(in)       :: ntiles
  real, dimension(:), intent(out) :: Lons, Lats
! local
  integer                  :: mark
  integer                  :: dum, n, nt
  character*128            :: dumstr
  real                     :: rdum1, rdum2
  
  
  open(10, file=tilefile, form="formatted",status='old')

  read(10, fmt=*) nt
  read(10, fmt=*) dum
  read(10, fmt=*) dumstr
  read(10, fmt=*) dum
  read(10, fmt=*) dum
  read(10, fmt=*) dumstr
  read(10, fmt=*) dum
  read(10, fmt=*) dum
  do
     read(10, fmt=*) mark
     if(mark .eq. 0) exit
  end do
  backspace(10)

  do n=1,ntiles
    read(10, fmt=*) mark, dum, rdum1, rdum2
    if(mark .eq. 0) then
        Lons(n) = rdum1
        Lats(n) = rdum2
    endif
  enddo
  !do n=1,10 
  !  print*,  Lons(n), Lats(n)
  !enddo

  close(10)

end subroutine GetTileLonLats

!=======================================================================
!BOP
!
! !IROUTINE: set_state_var - initialize single-category state variables
!
! !INTERFACE:
!
      subroutine set_state_var (nx_block, ny_block, &
                                icells,             &
                                indxi,       indxj, &
                                tile_lat,  tile_lon,         &
                                ice_con,  tsc,      &
                                aicen,    trcrn,    &
                                vicen,    vsnon,    &
                                eicen,    esnon) 
!
! !DESCRIPTION:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block , & ! block dimensions
         icells     

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice

      real (kind=real_kind), dimension (nx_block,ny_block), intent(in) :: &
         tile_lat, tile_lon, & 
         ice_con, tsc    


      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(out) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), &
         intent(out) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(out) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(out) :: &
         esnon     ! energy of melting for each ice layer  (J/m^2)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k           , & ! ice layer index
         n           , & ! thickness category index
         nc          , & ! thickness category index
         it              ! tracer index

      real (kind=dbl_kind) :: &
         slope, Ti, hi,  zn, hbar, &
         ainit(ncat), &
         hinit(ncat)

      real (kind=dbl_kind), dimension (nx_block,ny_block)  :: &
         ice_cov

      real (kind=dbl_kind), dimension(nilyr+1) :: &
         salin       , & ! salinity (ppt)   
         Tmlt            ! melting temp, -depressT * salinity
                         ! nilyr + 1 index is for bottom surface

      integer (kind=int_kind), parameter :: &
         nt_Tsfc  =  1, & ! ice/snow surface temperature
         nt_iage  =  2, & ! volume-weighted ice age
         nt_volpn =  3    ! melt pond volume

      real (kind=dbl_kind) :: &
         hin_max(0:ncat) ! category limits (m)

      real (kind=dbl_kind), parameter :: &
        pi   = 3.14159265358979323846_dbl_kind,&! pi
        c0   = 0.0_dbl_kind, &
        c1   = 1.0_dbl_kind, &
        c1p5 = 1.5_dbl_kind, &
        c2   = 2.0_dbl_kind, &
        c3   = 3.0_dbl_kind, &
        c4   = 4.0_dbl_kind, &
        c5   = 5.0_dbl_kind, &
        c6   = 6.0_dbl_kind, &
        c8   = 8.0_dbl_kind, &
        c9   = 9.0_dbl_kind, &
        c10  = 10.0_dbl_kind, &
        c12  = 12.0_dbl_kind, &
        c15  = 15.0_dbl_kind, &
        c16  = 16.0_dbl_kind, &
        c20  = 20.0_dbl_kind, &
        c25  = 25.0_dbl_kind, &
        c100 = 100.0_dbl_kind, &
        c180 = 180.0_dbl_kind, &
        c360 = 360.0_dbl_kind, &
        c365 = 365.0_dbl_kind, &
        c3600= 3600.0_dbl_kind, &
        c1000= 1000.0_dbl_kind, &
        p001 = 0.001_dbl_kind, &
        p01  = 0.01_dbl_kind, &
        p1   = 0.1_dbl_kind, &
        p2   = 0.2_dbl_kind, &
        p4   = 0.4_dbl_kind, &
        p5   = 0.5_dbl_kind, &
        p6   = 0.6_dbl_kind, &
        p05  = 0.05_dbl_kind, &
        p15  = 0.15_dbl_kind, &
        p25  = 0.25_dbl_kind, &
        p75  = 0.75_dbl_kind, &
        p166 = c1/c6, &
        p333 = c1/c3, &
        p666 = c2/c3, &
        p111 = c1/c9, &
        p055 = p111*p5, &
        p027 = p055*p5, &
        p222 = c2/c9, &
        eps04  = 1.0e-4_dbl_kind, &
        eps11  = 1.0e-11_dbl_kind, &
        eps13  = 1.0e-13_dbl_kind, &
        eps16  = 1.0e-16_dbl_kind, &
        puny   = eps11, &
        bignum = 1.0e+30_dbl_kind, &
        pih    = p5*pi, &
        pi2    = c2*pi



      real (kind=dbl_kind), parameter :: &
         tnh = 1.0_dbl_kind,           & 
         tsh = 0.75_dbl_kind,          &
         hsno_init = 0.20_dbl_kind   , & ! initial snow thickness (m)
         edge_init_nh =  70._dbl_kind, & ! initial ice edge, N.Hem. (deg) 
         edge_init_sh = -60._dbl_kind    ! initial ice edge, S.Hem. (deg)

      real (kind=dbl_kind), parameter :: &
         cp_ice    = 2106._dbl_kind   ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = 4218._dbl_kind   ,&! specific heat of ocn    (J/kg/K)
                                        ! freshwater value needed for enthalpy
         rhoi      = 917.0_dbl_kind   ,&! density of ice (kg/m^3)
         Lsub      = 2.835e6_dbl_kind ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501e6_dbl_kind ,&! latent heat, vaporization freshwater (J/kg)
         Lfresh    = Lsub-Lvap        ,&! latent heat of melting of fresh ice (J/kg)
         depressT  = 0.054_dbl_kind   ,&! Tf:brine salinity ratio (C/ppt)
         Tf        = -1.8_dbl_kind    ,& ! freezing temp.
         Tsmelt    = 0.0_dbl_kind     ,&! melting temperature, snow top surface (C)
         saltmax   = 3.2_dbl_kind     ,& ! max salinity at ice base (ppt)
         nsal      = 0.407_dbl_kind   ,&
         msal      = 0.573_dbl_kind   ,&
         min_salin = 0.1_dbl_kind     ! threshold for brine pocket treatment

      real (kind=dbl_kind) :: &
           cc1, cc2, cc3, & ! parameters for kcatbound = 0
           x1           , &
           rn           , & ! real(n)
           rncat        , & ! real(ncat)
           d1           , & ! parameters for kcatbound = 1 (m)
           d2
      real (kind=real_kind) :: lono, lato



      lono = -164.2552
      lato = 59.18839

      rncat = real(ncat, kind=dbl_kind)

      ! linear remapping itd category limits
      cc1 = c3/rncat
      cc2 = c15*cc1
      cc3 = c3

      hin_max(0) = c0     ! minimum ice thickness, m
      do n = 1, ncat
         x1 = real(n-1,kind=dbl_kind) / rncat
         hin_max(n) = hin_max(n-1) &
                    + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
      enddo

      do k = 1, nilyr
         zn = (real(k,kind=dbl_kind)-p5) /  &
               real(nilyr,kind=dbl_kind)
         salin(k)=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
         Tmlt(k) = -salin(k)*depressT
      enddo
      salin(nilyr+1) = saltmax
      Tmlt(nilyr+1) = -salin(nilyr+1)*depressT


      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         ice_cov(i,j) = ice_con(i,j)
         if (ice_cov(i,j) .lt. eps04) ice_cov(i,j) = c0
         if (ice_cov(i,j) .gt. c1)    ice_cov(i,j) = c1
      enddo

      ! Initialize state variables.
      ! If restarting, these values are overwritten.

      do n = 1, ncat
        do ij = 1, icells
           i = indxi(ij)
           j = indxj(ij)
           aicen(i,j,n) = c0
           vicen(i,j,n) = c0
           vsnon(i,j,n) = c0
           trcrn(i,j,nt_Tsfc,n) = Tf  ! surface temperature
        enddo
      enddo
      eicen(:,:,:) = c0
      esnon(:,:,:) = c0



      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)


         !--------------------------------------------------------------
         ! Place ice where ice concentration > .0001
         !--------------------------------------------------------------
         if (ice_cov(i,j) >= eps04) then

            hi = -1.0e10_dbl_kind
            !----------------------------------------------------------
            ! Set ice thickness in each hemisphere
            !----------------------------------------------------------
            if(tile_lat(i,j) > 30.0_dbl_kind) then
              hi  = tnh
            else if(tile_lat(i,j) < -30.0_dbl_kind) then
              hi  = tsh
            end if

            do nc = 1,ncat
              if(hin_max(nc-1) < hi .and. hi < hin_max(nc)) then

                  aicen(i,j,nc) = ice_cov(i,j)
                  vicen(i,j,nc) = hi*aicen(i,j,nc) 
                  trcrn(i,j,nt_Tsfc,nc) = min(Tsmelt, tsc(i,j) - Tffresh) !deg C
                  if(abs(tile_lon(i,j)-lono)<1.e-4 .and. abs(tile_lat(i,j)-lato)<1.e-4) then
                       print*,  tsc(i,j), trcrn(i,j,nt_Tsfc,nc)
                  endif 

                  do k = 1, nilyr

                     ! assume linear temp profile and compute enthalpy
                     slope = Tf - trcrn(i,j,nt_Tsfc,nc)
                     Ti = trcrn(i,j,nt_Tsfc,nc) &
                        + slope*(real(k,kind=dbl_kind)-p5) &
                                /real(nilyr,kind=dbl_kind)

                     eicen(i,j,ilyr1(nc)+k-1) = &
                          -(rhoi * (cp_ice*(Tmlt(k)-Ti) &
                          + Lfresh*(c1-Tmlt(k)/Ti) - cp_ocn*Tmlt(k))) &
                          * vicen(i,j,nc)/real(nilyr,kind=dbl_kind)
                     if(abs(tile_lon(i,j)-lono)<1.e-4 .and. &
                        abs(tile_lat(i,j)-lato)<1.e-4) then
                       print*, k, Ti
                     endif 
                  enddo               ! nilyr
              endif  
            enddo !ncat

         end if ! ice_cov(i,j) >= eps04
        enddo                  ! icells

      end subroutine set_state_var


end program cv_SaltRestart
