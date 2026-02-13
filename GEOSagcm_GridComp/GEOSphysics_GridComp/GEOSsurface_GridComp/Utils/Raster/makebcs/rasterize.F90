#include "MAPL_ErrLog.h"

module LogRectRasterizeMod

  use MAPL_SORTMOD
  use MAPL_ExceptionHandling  
  use MAPL_Constants,          only: PI=>MAPL_PI_R8
  use MAPL
  use catch_constants,        ONLY: CATCH_N_PFAFS
  use, intrinsic :: iso_fortran_env, only: INT32, REAL64
  implicit none
  private

!BOP

! !MODULE: LogRectRasterizeMod

! !DESCRIPTION: Container module for LRRasterize subroutine. It
!   has two public quantities: the subroutine itself and the value
!   of RASTERUNDEF

!EOP

  public LRRasterize
! public ReadRaster
  public WriteRaster
  public Writetiling
  public Sorttiling
  public MAPL_UNDEF_R8
  ! SRTM_maxcat = number of Pfafstetter catchments defined in raster file produced by Kristine Version in 2013
  !               (based on DEMs from 3.0-arcsec HydroSHEDS/SRTM south of 60N, 
  !                                   7.5-arcsec GMTED2010 north of 60N, and 
  !                                   CGIAR/SRTM where HydroSHEDS/SRTM is undefined [typically islands])  

  ! -------------------------------------------------------------------------------------------------------------
  INTEGER, PARAMETER, public:: SRTM_maxcat = CATCH_N_PFAFS

  integer,      parameter :: PUSHLEFT      = 10000
  real(REAL64), parameter :: Zero          = 0.0d0

  integer,      parameter :: NX            = 8640
  integer,      parameter :: NY            = 4320
  real(REAL64), parameter :: MAPL_UNDEF_R8 = 1.0D15

  real(REAL64)            :: garea_
  integer                 :: ctg_

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

! subroutine ReadRaster(File, Raster, Zip)
!  character*(*), intent(IN) :: File
!  integer,       intent(IN) :: Raster(:,:)
!  logical, optional         :: Zip
!
!  logical :: DoZip, Opened
!  integer :: nx, ny
!
!  nx = size(Raster,1)
!  ny = size(Raster,2)
!
!  if(present(Zip)) then
!     DoZip = Zip
!  else
!     DoZip = .false.
!  endif
!
!  if(DoZip) then
!     print *, "Reading zipped raster files not supported"
!     call exit(1)
!  else
!     call READRST(RASTER(1,1),nx,ny,trim(FILE)//CHAR(0))
!  end if
!
!  return
! end subroutine ReadRaster



subroutine SortTiling(Raster,rTable,iTable)
  integer, intent(INOUT) :: Raster(:,:), iTable(0:,:)
  real(REAL64),   intent(INOUT) :: rTable(:,:)

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

! --------------------------------------------------------------------------------------------

subroutine WriteTilingIR(File, GridName, im, jm, ipx, nx, ny, iTable, rTable, Zip, Verb, rc)

  ! Write ASCII tile file 
  !
  ! We have ascii tile files that support either 1 or 2 grids. The most typcal case for GEOS is a tile file with
  ! 2 grids (although to generate the final file, we generate a lot of intermediate tile files with a single grid).
  ! For the 2-grid tile file, the pfaf index number should be in column 9 (basically that file has 3 groups, 
  ! each of which has 4 columns: the "global" info part (type, area, lat, lon), and then for each grid we have 
  ! (index1 (i.e. "i"), index2 (i.e. "j"), weight, dummy). Here "dummy" is a variable, used internally for 
  ! bookkeeping purposes, but it is totally ignored by GEOS, MAPL, etc. So, for the typical case, ATM and OCN 
  ! grids, columns 1-4 represent the global variables, then the next 4 columns refer to the ATM grid (but this 
  ! is to a large extend an artifact of the ordering of the "combine" calls that generate the final tile file). 
  ! Then for type=0 (i.e., "ocean") the last 4 columns are the i, j, weight, dummy of the ocean grid. 
  ! But for type=100 (i.e., land) the convention is the first index, i.e. column 9, is the pfaf index
  ! (that is, the index of the Pfafstetter hydrological catchment). 
  ! I do not think we use the content of column 10 anywhere in the model.
  ! So my bottom line is the pfaf index should be in column 9. If it appears in column 8, it won't do any harm 
  ! to the atmosphere, but we cannot use it properly to do river routing inside the land model.
  ! (From https://github.com/GEOS-ESM/GEOSgcm_GridComp/pull/1028#issuecomment-2599275578, lightly edited.)


  character*(*),     intent(IN) :: File
  character*(*),     intent(IN) :: GridName(:)
  integer,           intent(IN) :: nx,ny
  integer,           intent(IN) :: iTable(0:,:)
  real(REAL64),      intent(IN) :: rTable(:,:)
  integer,           intent(IN) :: IM(:), JM(:), ipx(:)
  logical, optional, intent(IN) :: Zip
  logical, optional, intent(IN) :: Verb
  integer, optional, intent(out) :: rc

! Table variables
!
!  iTable(0)    :: Surface type
!  iTable(1)    :: tile count 
!  iTable(2)    :: I_1 I of 1st grid
!  iTable(3)    :: J_1
!  iTable(4)    :: I_2 I of 2nd grid *OR* for land tiles: index of Pfafstetter catchment (see comment above)
!  iTable(5)    :: J_2
!  iTable(6)    :: kk_1 (dummy variable for internal bookkeeping)
!  iTable(7)    :: kk_2 (dummy variable for internal bookkeeping)
!
!  rTable(1)    :: sum of lons
!  rTable(2)    :: sum of lats
!  rTable(3)    :: area
!  rTable(4)    :: of first grid box area
!  rTable(5)    :: of 2nd   grid box area

  logical        :: DoZip, Opened
  integer        :: j, unit, ng, ip, l, i, k, ix
  character*1000 :: Line
  integer        :: ii(size(GridName)), jj(size(GridName)), kk(size(GridName))
  real(REAL64)   :: fr(size(GridName))
  real(REAL64)   :: xc, yc, area
  real(REAL64)   :: garea, ctg(size(Gridname))
  real(REAL64)   :: sphere, error
  integer        :: status, tmp_in1, tmp_in2, ncat
  logical        :: file_exists

  ip = size(iTable,2)
  ng = size(GridName)

  _ASSERT(IP==size(rTable,2),'needs informative message')
  _ASSERT(NG==size(IM),      'needs informative message')
  _ASSERT(NG==size(JM),      'needs informative message')

  if(present(Zip)) then
     DoZip = Zip
  else
     DoZip = .false.
  endif

! Open unit and write header

  if(DoZip) then
     call ZipBuff(trim(File)//char(0),trim(Line)//char(0),0,1,1)
     write(LINE,'(4I10)') ip, SRTM_maxcat, nx, ny
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
     write(UNIT,'(4I10)') ip, SRTM_maxcat, nx, ny             
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
#define LINE_VARIABLES  iTable(0,k),area,xc,yc, (ii(l),jj(l),fr(l),kk(l),l=1,ng)   ! for *land* tiles, ii(2) = index of Pfafstetter catchment
 
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
        ! kk = "dummy" variable, used internally for bookkeeping purposes, ignored by GEOS, MAPL, etc
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
     sphere = 4.*PI
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

end subroutine WriteTilingIR

! ----------------------------------------------------------------------------------

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
     write(LINE,'(4I10)') ip, SRTM_maxcat, nx, ny
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
     write(UNIT,'(4I10)') ip,SRTM_maxcat, nx, ny             
     write(UNIT,'(I10)' ) 1
     ix=index(GridName   ,'/',Back=.true.)
     write(UNIT, *       ) trim(GridName   (ix+1:))
     write(UNIT,'(I10)' ) IM
     write(UNIT,'(I10)' ) JM
  end if

  return

end subroutine OpenTiling

! --------------------------------------------------------------------------

subroutine WriteLine(File, Unit, iTable, rTable, k, Zip, Verb)
  character*(*),     intent(IN) :: File
  integer,           intent(IN) :: Unit, k
  integer,           intent(IN) :: iTable(0:)
  real(REAL64),             intent(IN) :: rTable(:)
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
  real(REAL64)   :: fr
  real(REAL64)   :: xc, yc, area

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
  real(REAL64)   :: sphere, error
  character*1000 :: Line

  Line=""

  if(present(Zip)) then
     DoZip = Zip
  else
     DoZip = .false.
  endif

  if(present(Verb)) then
     sphere = 4.*PI
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


