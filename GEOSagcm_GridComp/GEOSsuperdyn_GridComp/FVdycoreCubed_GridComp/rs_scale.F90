      program  main
      use MAPL_IOMod
      use MAPL_ConstantsMod, only: MAPL_PSDRY
      use pFIO
      implicit none

! ************************************************************************
! ************************************************************************
! ****                                                                ****
! ****   Program to insert DRY MASS value of 983.05 mb into Restarts  ****
! ****   (to within 1e-10 Pa)                                         ****
! ****                                                                ****
! ************************************************************************
! ************************************************************************

      character*1024 :: Usage="rs_scale.x  fv_internal_rst  moist_internal_rst"
      character*1024 :: dynrst, mstrst
      character*1024 :: arg(2)

      integer headr1(6)
      integer headr2(5)
      integer nymd,nhms
      integer im,jm,lm
      integer i,j,L,n,rc,nargs,iargc
      integer iter

! restart variables
! -----------------
      real*8,   allocatable ::    u(:,:,:)
      real*8,   allocatable ::    v(:,:,:)
      real*8,   allocatable ::   th(:,:,:)
      real*8,   allocatable ::   pk(:,:,:)
      real*8,   allocatable ::   dp(:,:,:)
      real*8,   allocatable ::  pke(:,:,:)
      real*8,   allocatable ::  ple(:,:,:)
      real*8,   allocatable :: dum8(:,:)
      real*8,   allocatable ::   ak(:)
      real*8,   allocatable ::   bk(:)
      real*8,   allocatable :: delz(:,:,:)
      real*8,   allocatable ::    w(:,:,:)

      real,   allocatable ::   qv(:,:,:)
      real,   allocatable :: qlls(:,:,:)
      real,   allocatable :: qlcn(:,:,:)
      real,   allocatable :: cfls(:,:,:)
      real,   allocatable :: cfcn(:,:,:)
      real,   allocatable :: qils(:,:,:)
      real,   allocatable :: qicn(:,:,:)
      real*8, allocatable :: area(:,:)
      real,   allocatable :: dum4(:,:)

      real*8, allocatable ::    gsum(:,:)
      real*8, allocatable ::    qsum(:,:)
      real*8, allocatable ::   psold(:,:)
      real*8, allocatable ::   psnew(:,:)
      real*8, allocatable :: pdryold(:,:)
      real*8, allocatable :: pdrynew(:,:)

      real*8, parameter   ::    pdry_ave = MAPL_PSDRY
      real*8              :: pdryold_ave
      real*8              :: pdrynew_ave
      real*8              :: pdrydif_ave
      real*8, parameter   :: eps = epsilon(1.0d-10)
      real                :: kappa

      type(Netcdf4_FileFormatter) :: InDyn,OutDyn,InMoist,OutMoist
      type(FileMetadata) :: InCfgDyn,OutCfgDyn,InCfgMoist,OutCfgMoist
      integer         :: filetypeDyn, filetypeMoist, nVars,nVarsMoist

      kappa = 2.0/7.0

! **********************************************************************

  nargs = iargc()
  if( nargs<2 .or. nargs>3 ) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     call exit(1)
  end if

  do n=1,nargs
  call getarg(n,arg(n))
  enddo

! Open Dynamics and Moist Internal Restarts
! -----------------------------------------
     read(arg(1),'(a)') dynrst
     read(arg(2),'(a)') mstrst

     call MAPL_NCIOGetFileType(dynrst,filetypeDyn)
     call MAPL_NCIOGetFileType(mstrst,filetypeMoist)
     if (filetypeDyn /= filetypeMoist) stop

     if (filetypeDyn == 0) then

        call InDyn%open(dynrst,pFIO_READ,rc=rc)
        InCfgDyn=InDyn%read(rc=rc)
        im = InCfgDyn%get_dimension('lon',rc=rc)
        jm = InCfgDyn%get_dimension('lat',rc=rc)
        lm = InCfgDyn%get_dimension('lev',rc=rc)
        call MAPL_IOCountNonDimVars(InCfgDyn,nVars)
        call InMoist%open(mstrst,pFIO_READ,rc=rc)
        InCfgMoist = InMoist%read(rc=rc)
        call MAPL_IOCountNonDimVars(InCfgMoist,nVarsMoist)

        nymd=0
        nhms=0

     else

        open(unit=10, file=trim(dynrst), form='unformatted')
        open(unit=20, file=trim(mstrst), form='unformatted')
 
! **********************************************************************
! ****                  Read dycore internal Restart                ****
! **********************************************************************

        read (10) headr1
        read (10) headr2

        nymd = headr1(1)*10000 &
             + headr1(2)*100   &
             + headr1(3)
        nhms = headr1(4)*10000 &
             + headr1(5)*100   &
             + headr1(6)

        im = headr2(1)
        jm = headr2(2)
        lm = headr2(3)

     end if

      print *
      print *, '   dyn restart filename: ',trim(dynrst)
      print *, ' moist restart filename: ',trim(mstrst)
      print *, '             resolution: ',im,jm,lm
      print *, '                   date: ',nymd,nhms
      print *

      allocate (   ak(lm+1)       )
      allocate (   bk(lm+1)       )
      allocate (    u(im,jm,lm)   )
      allocate (    v(im,jm,lm)   )
      allocate (   th(im,jm,lm)   )
      allocate (   pk(im,jm,lm)   )
      allocate (  ple(im,jm,lm+1) )
      allocate (  pke(im,jm,lm+1) )

      if (filetypeDyn == 0) then

         call MAPL_VarRead(InDyn,"AK",ak)
         call MAPL_VarRead(InDyn,"BK",bk)
         do L=1,lm
            call MAPL_VarRead(InDyn,"U",u(:,:,L),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarRead(InDyn,"V",v(:,:,L),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarRead(InDyn,"PT",th(:,:,L),lev=l)
         enddo
         do L=1,lm+1
            call MAPL_VarRead(InDyn,"PE",ple(:,:,L),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarRead(InDyn,"PKZ",pk(:,:,L),lev=l)
         enddo
         ! check if we have delz and w
         if (nvars == 9) then
            allocate( delz(im,jm,lm ) )
            allocate(    w(im,jm,lm ) )
            do L=1,lm
               call MAPL_VarRead(InDyn,"DZ",delz(:,:,L),lev=l)
            enddo
            do L=1,lm
               call MAPL_VarRead(InDyn,"W",w(:,:,L),lev=l)
            enddo
         end if

      else

         read (10) ak
         read (10) bk

         do L=1,lm
            read(10)  u(:,:,L)
         enddo
         do L=1,lm
            read(10)  v(:,:,L)
         enddo
         do L=1,lm
            read(10)  th(:,:,L)
         enddo
         do L=1,lm+1
            read(10) ple(:,:,L)
         enddo
         do L=1,lm
            read(10)  pk(:,:,L)
         enddo

         close (10)

      endif

      allocate ( dp(im,jm,lm) )
      do L=1,lm
         dp(:,:,L) = ple(:,:,L+1)-ple(:,:,L)
      enddo

! **********************************************************************
! ****                   Read moist internal Restart                ****
! **********************************************************************

      allocate (   qv(im,jm,lm) )
      allocate ( qlls(im,jm,lm) )
      allocate ( qlcn(im,jm,lm) )
      allocate ( cfls(im,jm,lm) )
      allocate ( cfcn(im,jm,lm) )
      allocate ( qils(im,jm,lm) )
      allocate ( qicn(im,jm,lm) )

      if (filetypeMoist == 0) then

         do L=1,lm
            call MAPL_VarRead(InMoist,"Q",qv(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarRead(InMoist,"QLLS",qlls(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarRead(InMoist,"QLCN",qlcn(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarRead(InMoist,"CLLS",cfls(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarRead(InMoist,"CLCN",cfcn(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarRead(InMoist,"QILS",qils(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarRead(InMoist,"QICN",qicn(:,:,l),lev=l)
         enddo

      else

         do L=1,lm
            read(20)   qv(:,:,L)
         enddo
         do L=1,lm
            read(20) qlls(:,:,L)
         enddo
         do L=1,lm
            read(20) qlcn(:,:,L)
         enddo
         do L=1,lm
            read(20) cfls(:,:,L)
         enddo
         do L=1,lm
            read(20) cfcn(:,:,L)
         enddo
         do L=1,lm
            read(20) qils(:,:,L)
         enddo
         do L=1,lm
            read(20) qicn(:,:,L)
         enddo

      close (20)

      end if

! **********************************************************************
! ****                 Compute/Import Grid-Cell Area                ****
! **********************************************************************

           allocate ( area(im,jm) )
      call Get_Areas ( area,im,jm )

! **********************************************************************
! ****                       Constrain PDRY                         ****
! **********************************************************************

      allocate ( gsum(5,2)   )
      allocate ( qsum(im,jm) )

      allocate (   psold(im,jm) )
      allocate (   psnew(im,jm) )
      allocate ( pdryold(im,jm) )
      allocate ( pdrynew(im,jm) )

      pdrydif_ave = 1.0d0
             iter = 1

      do while ( dabs( pdrydif_ave ).gt.eps .and. iter.le.20 )

! --------------------------------------
      
      do n=1,5
         qsum = 0.0_8
         do L=1,lm
         if( n.eq.1 ) qsum = qsum +   qv(:,:,L)*dp(:,:,L)
         if( n.eq.2 ) qsum = qsum + qlls(:,:,L)*dp(:,:,L)
         if( n.eq.3 ) qsum = qsum + qlcn(:,:,L)*dp(:,:,L)
         if( n.eq.4 ) qsum = qsum + qils(:,:,L)*dp(:,:,L)
         if( n.eq.5 ) qsum = qsum + qicn(:,:,L)*dp(:,:,L)
         enddo
         call AreaMean( qsum, area, gsum(n,1), im,jm )
      enddo

      qsum = 0.0_8
      do L=1,lm
      qsum = qsum + (  qv(:,:,L) + &
                     qlls(:,:,L) + &
                     qlcn(:,:,L) + &
                     qils(:,:,L) + &
                     qicn(:,:,L) ) * dp(:,:,L)
      enddo
        psold = ple(:,:,lm+1)
      pdryold = ple(:,:,lm+1) - qsum  ! Subtract Total Water Content

      call AreaMean( pdryold, area, pdryold_ave, im,jm )

      pdrynew = pdryold * ( pdry_ave/pdryold_ave )

        psnew = pdrynew + qsum

      do L=1,lm+1
      do j=1,jm
      do i=1,im
       ple(i,j,L) = ple(i,j,L) + bk(L)*( psnew(i,j)-psold(i,j) )
      enddo
      enddo
      enddo

      pke = ple**kappa
      do L=1,lm
      dp(:,:,L) =             ple(:,:,L+1)-ple(:,:,L)
      pk(:,:,L) =            (pke(:,:,L+1)-pke(:,:,L)) &
                / ( kappa*log(ple(:,:,L+1)/ple(:,:,L)) )
      enddo

! --------------------------------------

      do n=1,5
         qsum = 0.0_8
         do L=1,lm
         if( n.eq.1 ) qsum = qsum +   qv(:,:,L)*dp(:,:,L)
         if( n.eq.2 ) qsum = qsum + qlls(:,:,L)*dp(:,:,L)
         if( n.eq.3 ) qsum = qsum + qlcn(:,:,L)*dp(:,:,L)
         if( n.eq.4 ) qsum = qsum + qils(:,:,L)*dp(:,:,L)
         if( n.eq.5 ) qsum = qsum + qicn(:,:,L)*dp(:,:,L)
         enddo
         call AreaMean( qsum, area, gsum(n,2), im,jm )
      enddo

        qv =   qv * ( gsum(1,1)/gsum(1,2) )
      qlls = qlls * ( gsum(2,1)/gsum(2,2) )
      qlcn = qlcn * ( gsum(3,1)/gsum(3,2) )
      qils = qils * ( gsum(4,1)/gsum(4,2) )
      qicn = qicn * ( gsum(5,1)/gsum(5,2) )

      qsum = 0.0_8
      do L=1,lm
      qsum = qsum + (  qv(:,:,L) + &
                     qlls(:,:,L) + &
                     qlcn(:,:,L) + &
                     qils(:,:,L) + &
                     qicn(:,:,L) ) * dp(:,:,L)
      enddo
      pdrynew = ple(:,:,lm+1) - qsum  ! Subtract Total Water Content

      call AreaMean( pdrynew, area, pdrynew_ave, im,jm )

      pdrydif_ave = pdrynew_ave - pdryold_ave

      write(6,1001) pdrynew_ave/100,pdryold_ave/100,pdrynew_ave/pdryold_ave,pdrydif_ave/100
 1001 format(1x,'PSDRY_NEW: ',g21.14,'  PSDRY_OLD: ',g21.14,'  RATIO: ',g25.18,'  DIF: ',g21.14)

      iter = iter + 1
      enddo

! **********************************************************************
! ****                 Write dycore internal Restart                ****
! **********************************************************************

      if (filetypeDyn == 0) then

         dynrst = trim(dynrst) // '.scaled'
         OutCfgDyn=InCfgDyn
         call OutDyn%create(dynrst,rc=rc)
         call OutDyn%write(OutCfgDyn,rc=rc)

         call MAPL_VarWrite(OutDyn,"AK",ak)
         call MAPL_VarWrite(OutDyn,"BK",bk)
         do L=1,lm
            call MAPL_VarWrite(OutDyn,"U",u(:,:,L),lev=L)
         enddo
         do L=1,lm
            call MAPL_VarWrite(OutDyn,"V",v(:,:,L),lev=L)
         enddo
         do L=1,lm
            call MAPL_VarWrite(OutDyn,"PT",th(:,:,L),lev=L)
         enddo
         do L=1,lm+1
            call MAPL_VarWrite(OutDyn,"PE",ple(:,:,L),lev=L)
         enddo
         do L=1,lm
            call MAPL_VarWrite(OutDyn,"PKZ",pk(:,:,L),lev=L)
         enddo
         if (nvars==9) then
            do L=1,lm
               call MAPL_VarWrite(OutDyn,"DZ",delz(:,:,L),lev=L)
            enddo
            do L=1,lm
               call MAPL_VarWrite(OutDyn,"W",w(:,:,L),lev=L)
            enddo
         end if

      else

         allocate ( dum8(im,jm) )

         open (10,file=trim(dynrst),form='unformatted',access='sequential')

         dynrst = trim(dynrst) // '.scaled'
         print *
         print *, 'Creating GEOS-5 fvcore_internal_restart: ',trim(dynrst)
         open (20,file=trim(dynrst),form='unformatted',access='sequential')

         read (10) headr1
         read (10) headr2
         read (10) ak
         read (10) bk

         write(20) headr1
         write(20) headr2
         write(20) ak
         write(20) bk

             do L=1,lm
                read (10) dum8
                write(20) u(:,:,L)
             enddo
             do L=1,lm
                read (10) dum8
                write(20) v(:,:,L)
             enddo
             do L=1,lm
                read (10) dum8
                write(20) th(:,:,L)
             enddo
             do L=1,lm+1
                read (10) dum8
                write(20) ple(:,:,L)
             enddo
             do L=1,lm
                read (10) dum8
                write(20) pk(:,:,L)
             enddo

                   rc =  0
         do while (rc.eq.0)
            read (10,iostat=rc)     dum8
            if( rc.eq.0 ) write(20) dum8
         enddo

         close (10)
         close (20)

      end if

! **********************************************************************
! ****                  Write moist internal Restart                ****
! **********************************************************************

      if (fileTypeMoist == 0) then

         mstrst = trim(mstrst) // '.scaled'

         OutCfgMoist=InCfgMoist
         call OutMoist%create(mstrst,rc=rc)
         call OutMoist%write(OutCfgMoist,rc=rc)

         do L=1,lm
            call MAPL_VarWrite(OutMoist,"Q",qv(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarWrite(OutMoist,"QLLS",qlls(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarWrite(OutMoist,"QLCN",qlcn(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarWrite(OutMoist,"CLLS",cfls(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarWrite(OutMoist,"CLCN",cfcn(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarWrite(OutMoist,"QILS",qils(:,:,l),lev=l)
         enddo
         do L=1,lm
            call MAPL_VarWrite(OutMoist,"QICN",qicn(:,:,l),lev=l)
         enddo
         ! check if we have ncpl and ncpi
         if (nVarsMoist == 9) then
            allocate( dum4(im,jm) )
               
            do L=1,lm
               call MAPL_VarRead(InMoist,"NCPL",dum4,lev=l)
               call MAPL_VarWrite(OutMoist,"NCPL",dum4,lev=l)
            enddo
            do L=1,lm
               call MAPL_VarRead(InMoist,"NCPI",dum4,lev=l)
               call MAPL_VarWrite(OutMoist,"NCPI",dum4,lev=l)
            enddo

         end if

      else

         allocate ( dum4(im,jm) )

         open (10,file=trim(mstrst),form='unformatted',access='sequential')

         mstrst = trim(mstrst) // '.scaled'
         print *
         print *, 'Creating GEOS-5  moist_internal_restart: ',trim(mstrst)
         open (20,file=trim(mstrst),form='unformatted',access='sequential')

         do L=1,lm
            read (10) dum4
            write(20)   qv(:,:,L)
         enddo
         do L=1,lm
            read (10) dum4
            write(20) qlls(:,:,L)
         enddo
         do L=1,lm
            read (10) dum4
            write(20) qlcn(:,:,L)
         enddo
         do L=1,lm
            read (10) dum4
            write(20) cfls(:,:,L)
         enddo
         do L=1,lm
            read (10) dum4
            write(20) cfcn(:,:,L)
         enddo
         do L=1,lm
            read (10) dum4
            write(20) qils(:,:,L)
         enddo
         do L=1,lm
            read (10) dum4
            write(20) qicn(:,:,L)
         enddo

                   rc =  0
         do while (rc.eq.0)
            read (10,iostat=rc)     dum4
            if( rc.eq.0 ) write(20) dum4
         enddo

         close (10)
         close (20)

      end if

      stop
      end

! --------------------------------------------------------------

      subroutine Get_Areas ( area,im,jm )
      use ESMF
      use constants_mod, only: cnst_radius=>radius
      use fv_grid_utils_mod, only: get_area
      use fv_arrays_mod,     only: R_GRID
      implicit   none
      integer    im,jm
      real*8  area(im,jm)

      real*8, allocatable :: cosp(:)
      real*8, allocatable ::   da(:,:)
      real*8              :: pi,dl,dp,acap,qdum
      integer i,j
      real(ESMF_KIND_R8), allocatable :: grid_lons(:,:,:)
      real(ESMF_KIND_R8), allocatable :: grid_lats(:,:,:)
      real(kind=R_GRID) :: p1(2),p2(2),p3(2),p4(2)
      integer :: jj,k

      real(kind=R_GRID), parameter :: radius = cnst_radius

#define R8  8
interface
  subroutine AppCSEdgeCreateF(IM_World, LonEdge, LatEdge, LonCenter, LatCenter, rc)
    integer,            intent(in   ) :: IM_World
    real(R8), intent(inout) :: LonEdge(IM_World+1,IM_World+1,6)
    real(R8), intent(inout) :: LatEdge(IM_World+1,IM_World+1,6)
    real(R8), intent(inout), optional :: LonCenter(IM_World,IM_World,6)
    real(R8), intent(inout), optional :: LatCenter(IM_World,IM_World,6)
    integer, optional,  intent(out  ) :: rc
  end subroutine AppCSEdgeCreateF
end interface

      ! LatLon Area
      ! -----------
      if( jm.ne.im*6 ) then
         allocate( da(im,jm) )
         allocate(  cosp(jm) )

         pi = 4.0D0 * Datan(1.0D0)
         dl = 2.0D0 * pi/ im
         dp =         pi/(jm-1)

         do j=2,jm-1
            cosp(j) = Dcos( -pi/2 + (j-1)*dp )
         enddo

         do j=1,jm
            if ( j == 1  ) then
                 dA(:,j) = 2*pi*(1-Dcos(dp/2))/im
            else if ( j == jm ) then
                 dA(:,j) = 2*pi*(1-Dcos(dp/2))/im
            else
                 dA(:,j) = cosp(j) * dl*dp
            endif
         enddo

         area = dA

         deallocate( da   )
         deallocate( cosp )

      else
      ! Cube Area
      ! ---------
        allocate(grid_lons(im+1,im+1,6)) 
        allocate(grid_lats(im+1,im+1,6))
        call AppCSEdgeCreateF(im,grid_lons,grid_lats)
        do k=1,6
           do j=1,im
              do i=1,im
                 jj=j + (k-1)*im
                 p1(1)=grid_lons(i,j,k)
                 p1(2)=grid_lats(i,j,k)
                 p2(1)=grid_lons(i,j+1,k)
                 p2(2)=grid_lats(i,j+1,k)
                 p3(1)=grid_lons(i+1,j,k)
                 p3(2)=grid_lats(i+1,j,k)
                 p4(1)=grid_lons(i+1,j+1,k)
                 p4(2)=grid_lats(i+1,j+1,k)
                 area(i,jj)=get_area(p1,p2,p3,p4,radius)
              enddo
           enddo
        enddo
      endif

      return
      end subroutine Get_Areas

! --------------------------------------------------------------

      subroutine AreaMean ( q,area,qave,im,jm )
      implicit     none
      integer      im,jm
      real*8  area(im,jm)
      real*8     q(im,jm)
      real*8     qave
      real*8     qdum

      integer i,j

      qave = 0.0_8
      qdum = 0.0_8
      do j=1,jm
         do i=1,im
            qave = qave + q(i,j)*area(i,j)
            qdum = qdum +        area(i,j)
         enddo
      enddo

      qave = qave / qdum

      return
      end subroutine AreaMean
