
#include "Raster.h"

module LogRectRasterizeMod

  use MAPL_SORTMOD

  implicit none
  private

!BOP

! !MODULE: LogRectRasterizeMod

! !DESCRIPTION: Container module for LRRasterize subroutine. It
!   has two public quantities: the subroutine itself and the value
!   of RASTERUNDEF

!EOP

  public LRRasterize
  public ReadRaster
  public WriteRaster
  public Writetiling
  public Sorttiling
  public Opentiling
  public Closetiling
  public WriteLine

  integer, parameter :: PUSHLEFT     = 10000
  REAL_  , parameter :: Zero         = 0.0
  REAL_  , parameter :: PI           = RASTER_PI

  integer, parameter :: NX           = 8640
  integer, parameter :: NY           = 4320


  REAL_   :: garea_
  integer :: ctg_
  
  interface LRRasterize
     module procedure LRRasterize2File
     module procedure LRRasterize2File0
  end interface

  interface WriteTiling
     module procedure WriteTilingIR
  end interface

contains

#include "rasterize.H"

#define MESH
#include "rasterize.H"
#undef MESH


subroutine WriteRaster(File, Raster, Zip)
  character*(*), intent(IN) :: File
  integer,       intent(IN) :: Raster(:,:)
  logical, optional         :: Zip

  integer :: nx, ny, DoZip

  nx = size(Raster,1)
  ny = size(Raster,2)

  DoZip = 0.

  if(present(Zip)) then
     if(Zip) DoZip =1
  endif

  call WRITERST(RASTER(1,1),nx,ny,trim(FILE)//CHAR(0),DoZip)

  return
end subroutine WriteRaster





subroutine ReadRaster(File, Raster, Zip)
  character*(*), intent(IN) :: File
  integer,       intent(IN) :: Raster(:,:)
  logical, optional         :: Zip

  logical :: DoZip, Opened
  integer :: nx, ny

  nx = size(Raster,1)
  ny = size(Raster,2)

  if(present(Zip)) then
     DoZip = Zip
  else
     DoZip = .false.
  endif

  if(DoZip) then
     print *, "Reading zipped raster files not supported"
     call exit(1)
  else
     call READRST(RASTER(1,1),nx,ny,trim(FILE)//CHAR(0))
  end if

  return
end subroutine ReadRaster



subroutine SortTiling(Raster,rTable,iTable)
  integer, intent(INOUT) :: Raster(:,:), iTable(0:,:)
  REAL_,   intent(INOUT) :: rTable(:,:)

  integer,   dimension(size(iTable,2)) :: old, new
  integer*8, dimension(size(iTable,2)) :: key, key0
  integer :: ip, i, j, k

  ip = size(rTable,2)

  do k=1,ip
     old(k) = k
     new(k) = k
  end do

! Sort table so non-ocean is first, in ascending Pfafstetter order
! and ocean is in ascending latitude.

  
  if(size(iTable,1)==4) then
     key0 = iTable(2,:ip)
     where(iTable(0,:ip)==0) key0 = iTable(3,:ip)*10000 + key0 + 10000000
  else
     key0 = iTable(4,:ip)
     where(iTable(0,:ip)==0)
        key0 = iTable(5,:ip)*10000 + key0 + 10000000
     end where
     key0 = key0*100000000 + iTable(6,:ip) ! add atmos index
  end if

  key = key0
  call MAPL_Sort(key,old)

  key = key0
  call MAPL_Sort(key,rTable(:,:ip))

  key = key0
  call MAPL_Sort(key,iTable(:,:ip))

  key = old
  call MAPL_Sort(key,new)


  do j=1,size(raster,2)
     do i=1,size(raster,1)
        raster(i,j) = new(raster(i,j))
     end do
  end do

  return
end subroutine SortTiling

subroutine WriteTilingIR(File, GridName, im, jm, ipx, nx, ny, iTable, rTable, Zip, Verb)
  character*(*),     intent(IN) :: File
  character*(*),     intent(IN) :: GridName(:)
  integer,           intent(IN) :: nx,ny
  integer,           intent(IN) :: iTable(0:,:)
  REAL_,             intent(IN) :: rTable(:,:)
  integer,           intent(IN) :: IM(:), JM(:), ipx(:)
  logical, optional, intent(IN) :: Zip
  logical, optional, intent(IN) :: Verb

! Table variables
!
!  iTable(0)    :: Surface type
!  iTable(1)    :: tile count 
!  iTable(2)    :: I_1 I of first grid
!  iTable(3)    :: J_1
!  iTable(4)    :: I_2 I of 2nd   grid
!  iTable(5)    :: J_2
!
!  rTable(1)    :: sum of lons
!  rTable(2)    :: sum of lats
!  rTable(3)    :: area
!  rTable(4)    :: of first grid box area
!  rTable(5)    :: of 2nd   grid box area

  logical :: DoZip, Opened
  integer :: j, unit, ng, ip, l, i, k, ix
  character*1000 :: Line
  integer :: ii(size(GridName)), jj(size(GridName)), kk(size(GridName))
  REAL_   :: fr(size(GridName))
  REAL_   :: xc, yc, area
  REAL_   :: garea, ctg(size(Gridname))
  REAL_   :: sphere, error

  ip = size(iTable,2)
  ng = size(GridName)

  ASSERT_(IP==size(rTable,2))
  ASSERT_(NG==size(IM))
  ASSERT_(NG==size(JM))

  if(present(Zip)) then
     DoZip = Zip
  else
     DoZip = .false.
  endif

! Open unit and write header

  if(DoZip) then
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),0,1,1)
     write(LINE,'(3I10)') ip, nx, ny
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
     write(LINE,'(I10)' ) ng
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
     do i=1,ng
        ix=index(GridName(i),'/',Back=.true.)
        write(LINE, *      ) trim(GridName(i)(ix+1:))
        call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
        write(LINE,'(I10)' ) IM(i)
        call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
        write(LINE,'(I10)' ) JM(i)            
        call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
     end do
  else
     Unit   = 9
     Opened = .true.
     do while(Opened)
        Unit = Unit + 1
        Inquire(Unit=Unit,Opened=Opened)
     end do
     open (UNIT,file=trim(File), form='formatted', status='unknown')
     write(UNIT,'(3I10)') ip, nx, ny             
     write(UNIT,'(I10)' ) ng              
     do l=1,ng
        ix=index(GridName(l),'/',Back=.true.)
        write(UNIT, *       ) trim(GridName(l)(ix+1:))
        write(UNIT,'(I10)' ) IM(l)
        write(UNIT,'(I10)' ) JM(l)
     end do
  end if

! Write tile info, one line per tile.

#define LINE_FORMAT     '(I10,3E20.12,9(2I10,E20.12,I10))'
#define LINE_VARIABLES  iTable(0,k),area,xc,yc, (ii(l),jj(l),fr(l),kk(l),l=1,ng)

  garea = 0.0
  ctg   = 0.0

  do k=1,ip
     xc    = rTable(1,k)
     yc    = rTable(2,k)
     area  = rTable(3,k)

     garea = garea + area

     do l=0,ng-1
        ii(l+1) = iTable(2 +L*2,K)
        jj(l+1) = iTable(3 +L*2,K)
        if(ng==1) then
           kk(l+1) = K
        else
           kk(l+1) = iTable(6 +L,K)
        end if
        if(rTable(4+L,K)/=0.0) then
           fr (l+1) = area / rTable(4+L,K)
           ctg(l+1) = ctg(l+1) + fr(l+1)
        elseif(Verb) then
           print *, 'Grid ',l+1,' appears to be undefined for tile ',k
        endif
     end do

     if(DoZip) then
        write(Line   ,LINE_FORMAT) LINE_VARIABLES
        call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
     else
        write(UNIT,LINE_FORMAT) LINE_VARIABLES
     end if
  end do

  if(present(Verb)) then
     sphere = 4.*pi
     error  = (sphere-garea)/garea
     if(Verb) then
        print '(A,3e20.13)','Stats for the globe:',garea, sphere, error
        do l=1,ng
           print *, 'Grid ', L, ' has ', IPX(l),' tiles; ', &
                    'overlay accounts for ',ctg(L),' of them.'
        end do
     end if
  end if

! Close the file

  if(DoZip) then
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),0,1,1)
  else
     close(UNIT)
  end if

  return
end subroutine WriteTilingIR


subroutine OpenTiling(Unit, File, GridName, im, jm, ip, nx, ny, Zip, Verb)
  integer,           intent(OUT) :: Unit
  character*(*),     intent(IN) :: File
  character*(*),     intent(IN) :: GridName
  integer,           intent(IN) :: nx,ny
  integer,           intent(IN) :: IM, JM, ip
  logical, optional, intent(IN) :: Zip
  logical, optional, intent(IN) :: Verb

! Table variables
!
!  iTable(0)    :: Surface type
!  iTable(1)    :: tile count 
!  iTable(2)    :: I_1 I of first grid
!  iTable(3)    :: J_1
!  iTable(4)    :: I_2 I of 2nd   grid
!  iTable(5)    :: J_2
!
!  rTable(1)    :: sum of lons
!  rTable(2)    :: sum of lats
!  rTable(3)    :: area
!  rTable(4)    :: of first grid box area
!  rTable(5)    :: of 2nd   grid box area

  logical :: DoZip, Opened
  integer :: j, l, i, k, ix
  character*1000 :: Line


  if(present(Zip)) then
     DoZip = Zip
  else
     DoZip = .false.
  endif

  garea_ = 0.0
  ctg_   = 0.0

! Open unit and write header

  if(DoZip) then
     Unit = -1
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),0,1,1)
     write(LINE,'(3I10)') ip, nx, ny
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
     write(LINE,'(I10)' ) 1
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
     ix=index(GridName   ,'/',Back=.true.)
     write(LINE, *      ) trim(GridName   (ix+1:))
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
     write(LINE,'(I10)' ) IM
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
     write(LINE,'(I10)' ) JM            
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
  else
     Unit   = 9
     Opened = .true.
     do while(Opened)
        Unit = Unit + 1
        Inquire(Unit=Unit,Opened=Opened)
     end do
     open (UNIT,file=trim(File), form='formatted', status='unknown')
     write(UNIT,'(3I10)') ip, nx, ny             
     write(UNIT,'(I10)' ) 1
     ix=index(GridName   ,'/',Back=.true.)
     write(UNIT, *       ) trim(GridName   (ix+1:))
     write(UNIT,'(I10)' ) IM
     write(UNIT,'(I10)' ) JM
  end if

  return
end subroutine OpenTiling




subroutine WriteLine(File, Unit, iTable, rTable, k, Zip, Verb)
  character*(*),     intent(IN) :: File
  integer,           intent(IN) :: Unit, k
  integer,           intent(IN) :: iTable(0:)
  REAL_,             intent(IN) :: rTable(:)
  logical, optional, intent(IN) :: Zip
  logical, optional, intent(IN) :: Verb

! Table variables
!
!  iTable(0)    :: Surface type
!  iTable(1)    :: tile count 
!  iTable(2)    :: I_1 I of first grid
!  iTable(3)    :: J_1
!  iTable(4)    :: I_2 I of 2nd   grid
!  iTable(5)    :: J_2
!
!  rTable(1)    :: sum of lons
!  rTable(2)    :: sum of lats
!  rTable(3)    :: area
!  rTable(4)    :: of first grid box area
!  rTable(5)    :: of 2nd   grid box area

  logical :: DoZip
  character*1000 :: Line
  integer :: ii, jj
  REAL_   :: fr
  REAL_   :: xc, yc, area

  if(present(Zip)) then
     DoZip = Zip
  else
     DoZip = .false.
  endif

! Write tile info, one line per tile.


#undef LINE_FORMAT
#undef LINE_VARIABLES
#define LINE_FORMAT     '(I10,3E20.12,9(2I10,E20.12,I10))'
#define LINE_VARIABLES  iTable(0),area,xc,yc, ii,jj,fr,k


  xc    = rTable(1)
  yc    = rTable(2)
  area  = rTable(3)

  garea_ = garea_ + area

  ii = iTable(2)
  jj = iTable(3)

  if(rTable(4)/=0.0) then
     fr  = area / rTable(4)
     ctg_ = ctg_ + fr
  elseif(Verb) then
     print *, 'Grid  appears to be undefined at this tile'
  endif

  if(DoZip) then
     write(Line,LINE_FORMAT) LINE_VARIABLES
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),1,1,1)
  else
     write(UNIT,LINE_FORMAT) LINE_VARIABLES
  end if

  return
end subroutine WriteLine


subroutine CloseTiling(FIle, Unit, ip,  Zip, Verb)
  character*(*),     intent(IN) :: File
  integer,           intent(IN) :: Unit, ip
  logical, optional, intent(IN) :: Zip
  logical, optional, intent(IN) :: Verb

! Table variables
!
!  iTable(0)    :: Surface type
!  iTable(1)    :: tile count 
!  iTable(2)    :: I_1 I of first grid
!  iTable(3)    :: J_1
!  iTable(4)    :: I_2 I of 2nd   grid
!  iTable(5)    :: J_2
!
!  rTable(1)    :: sum of lons
!  rTable(2)    :: sum of lats
!  rTable(3)    :: area
!  rTable(4)    :: of first grid box area
!  rTable(5)    :: of 2nd   grid box area

  logical :: DoZip
  REAL_   :: sphere, error
  character*1000 :: Line

  Line=""

  if(present(Zip)) then
     DoZip = Zip
  else
     DoZip = .false.
  endif

  if(present(Verb)) then
     sphere = 4.*pi
     error  = (sphere-garea_)/garea_
     if(Verb) then
        print '(A,3e20.13)','Stats for the globe:',garea_, sphere, error
        print *, 'Grid  has ', IP,' tiles; ', &
                    'overlay accounts for ',ctg_,' of them.'
     end if
  end if

! Close the file

  if(DoZip) then
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),0,1,1)
  else
     close(UNIT)
  end if

  return
end subroutine CloseTiling





end module LogRectRasterizeMod


