#include "MAPL_ErrLog.h"

module LogRectRasterizeMod

  use MAPL_SORTMOD
  use MAPL_ExceptionHandling  
  use MAPL_Constants,          only: PI=>MAPL_PI_R8
  use MAPL
  use, intrinsic :: iso_fortran_env, only: INT32
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
  public WritetilingNC4
  public ReadTilingNC4
  public Sorttiling
  public Opentiling
  public Closetiling
  public WriteLine

  integer,      parameter :: PUSHLEFT      = 10000
  real(kind=8), parameter :: Zero          = 0.0

  integer,      parameter :: NX            = 8640
  integer,      parameter :: NY            = 4320
  real(kind=8), parameter :: MAPL_UNDEF_R8 = MAPL_UNDEF

  real(kind=8)            :: garea_
  integer                 :: ctg_
  INTEGER, PARAMETER, public:: SRTM_maxcat = 291284 
 
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
  real(kind=8),   intent(INOUT) :: rTable(:,:)

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

subroutine WriteTilingIR(File, GridName, im, jm, ipx, nx, ny, iTable, rTable, Zip, Verb, rc)
  character*(*),     intent(IN) :: File
  character*(*),     intent(IN) :: GridName(:)
  integer,           intent(IN) :: nx,ny
  integer,           intent(IN) :: iTable(0:,:)
  real(kind=8),      intent(IN) :: rTable(:,:)
  integer,           intent(IN) :: IM(:), JM(:), ipx(:)
  logical, optional, intent(IN) :: Zip
  logical, optional, intent(IN) :: Verb
  integer, optional, intent(out) :: rc

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

  logical        :: DoZip, Opened
  integer        :: j, unit, ng, ip, l, i, k, ix
  character*1000 :: Line
  integer        :: ii(size(GridName)), jj(size(GridName)), kk(size(GridName))
  real(kind=8)   :: fr(size(GridName))
  real(kind=8)   :: xc, yc, area
  real(kind=8)   :: garea, ctg(size(Gridname))
  real(kind=8)   :: sphere, error
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

end subroutine WriteTilingIR

subroutine WriteTilingNC4(File, GridName, im, jm, nx, ny, iTable, rTable, srtm, rc)
  character(*),      intent(IN) :: File
  character(*),      intent(IN) :: GridName(:)
  integer,           intent(IN) :: IM(:), JM(:)
  integer,           intent(IN) :: nx,ny
  integer,           intent(IN) :: iTable(:,0:)
  real(kind=8),      intent(IN) :: rTable(:,:)
  integer, optional, intent(in) :: srtm
  integer, optional, intent(out):: rc

  integer                       :: k, ll, ng, ip, status, srtm_maxcat_

  character(len=:), allocatable :: attr
  type (Variable)               :: v
  type (NetCDF4_FileFormatter)  :: formatter
  character(len=1)              :: str_num
  type (FileMetadata)           :: metadata
  integer,          allocatable :: II(:), JJ(:), KK(:)
  real(kind=8),     allocatable :: fr(:)
  logical                       :: EASE

  ng  = size(GridName)
  ip  = size(iTable,1)

  EASE = .false.
  if (index(GridName(1), 'EASE') /=0) EASE = .true.

  srtm_maxcat_ = SRTM_maxcat
  if (present(srtm)) srtm_maxcat_ = srtm

  call metadata%add_dimension('tile', ip)

  do ll = 1, ng
    if (ng == 1) then
      str_num = ''
    else
      write(str_num, '(i0)') ll
    endif

    attr = 'Grid'//trim(str_num)//'_Name'
    call metadata%add_attribute( attr, trim(GridName(ll)))
    attr = 'IM'//trim(str_num)
    call metadata%add_attribute( attr, IM(ll))
    attr = 'JM'//trim(str_num)
    call metadata%add_attribute( attr, JM(ll))
  enddo

  attr = 'raster_nx'
  call metadata%add_attribute( attr, nx)
  attr = 'raster_ny'
  call metadata%add_attribute( attr, ny)
  attr = 'SRTM_maxcat'
  call metadata%add_attribute( attr, srtm_maxcat_)

  v = Variable(type=PFIO_INT32, dimensions='tile')
  call v%add_attribute('units', '1')
  call v%add_attribute('long_name', 'tile_type')
  call metadata%add_variable('typ', v)

  v = Variable(type=PFIO_REAL64, dimensions='tile')
  call v%add_attribute('units', 'km2')
  call v%add_attribute('long_name', 'tile_area')
  call v%add_attribute("missing_value", MAPL_UNDEF_R8)
  call v%add_attribute("_FillValue", MAPL_UNDEF_R8)
  call metadata%add_variable('area', v)

  v = Variable(type=PFIO_REAL64, dimensions='tile')
  call v%add_attribute('units', 'degree')
  call v%add_attribute('long_name', 'tile_center_of_mass_longitude')
  call v%add_attribute("missing_value", MAPL_UNDEF_R8)
  call v%add_attribute("_FillValue", MAPL_UNDEF_R8)
  call metadata%add_variable('com_lon', v)
  
  v = Variable(type=PFIO_REAL64, dimensions='tile')
  call v%add_attribute('units', 'degree')
  call v%add_attribute('long_name', 'tile_center_of_mass_latitude')
  call v%add_attribute("missing_value", MAPL_UNDEF_R8)
  call v%add_attribute("_FillValue", MAPL_UNDEF_R8)
  call metadata%add_variable('com_lat', v)
 
  do ll = 1, ng
     if (ng == 1) then
        str_num = ''
     else
        write(str_num, '(i0)') ll
     endif

     v = Variable(type=PFIO_INT32, dimensions='tile')
     call v%add_attribute('units', '1')
     call v%add_attribute('long_name', 'GRID'//trim(str_num)//'_i_index_of_tile_in_global_grid')
     call metadata%add_variable('i_indg'//trim(str_num), v)

     v = Variable(type=PFIO_INT32, dimensions='tile')
     call v%add_attribute('units', '1')
     call v%add_attribute('long_name', 'GRID'//trim(str_num)//'_j_index_of_tile_in_global_grid')
     call metadata%add_variable('j_indg'//trim(str_num), v)

     v = Variable(type=PFIO_REAL64, dimensions='tile')
     call v%add_attribute('units', '1')
     call v%add_attribute('long_name', 'GRID'//trim(str_num)//'_area_fraction_of_tile_in_grid_cell')
     call v%add_attribute("missing_value", MAPL_UNDEF_R8)
     call v%add_attribute("_FillValue",    MAPL_UNDEF_R8)
     call metadata%add_variable('frac_cell'//trim(str_num), v)

     v = Variable(type=PFIO_INT32, dimensions='tile')
     call v%add_attribute('units', '1')
     call v%add_attribute('long_name', 'Pfafstetter_index')
     call metadata%add_variable('pfaf_index'//trim(str_num), v)
  enddo

  if ( .not. EASE ) then
     v = Variable(type=PFIO_REAL64, dimensions='tile')
     call v%add_attribute('units', '1')
     call v%add_attribute('long_name', 'fraction_of_Pfafstetter')
     call metadata%add_variable('frac_pfaf', v)
  endif

  v = Variable(type=PFIO_REAL64, dimensions='tile')
  call v%add_attribute('units', 'degree')
  call v%add_attribute('long_name', 'tile_minimum_longitude')
  call v%add_attribute("missing_value", MAPL_UNDEF_R8)
  call metadata%add_variable('min_lon', v)

  v = Variable(type=PFIO_REAL64, dimensions='tile')
  call v%add_attribute('units', 'degree')
  call v%add_attribute('long_name', 'tile_maximum_longitude')
  call v%add_attribute("missing_value", MAPL_UNDEF_R8)
  call metadata%add_variable('max_lon', v)

  v = Variable(type=PFIO_REAL64, dimensions='tile')
  call v%add_attribute('units', 'degree')
  call v%add_attribute('long_name', 'tile_minimum_latitude')
  call v%add_attribute("missing_value", MAPL_UNDEF_R8)
  call metadata%add_variable('min_lat', v)

  v = Variable(type=PFIO_REAL64, dimensions='tile')
  call v%add_attribute('units', 'degree')
  call v%add_attribute('long_name', 'tile_maximum_latitude')
  call v%add_attribute("missing_value", MAPL_UNDEF_R8)
  call metadata%add_variable('max_lat', v)

  v = Variable(type=PFIO_REAL64, dimensions='tile')
  call v%add_attribute('units', 'm')
  call v%add_attribute('long_name', 'tile_mean_elevation')
  call v%add_attribute("missing_value", MAPL_UNDEF_R8)
  call metadata%add_variable('elev', v)

  call formatter%create(File, mode=PFIO_NOCLOBBER, rc=status)
  call formatter%write(metadata,                   rc=status)
  call formatter%put_var('typ',     iTable(:,0),   rc=status)
  call formatter%put_var('area',    rTable(:,3),   rc=status)
  call formatter%put_var('com_lon', rTable(:,1),   rc=status)
  call formatter%put_var('com_lat', rTable(:,2),   rc=status)

  allocate(fr(ip))
  fr = MAPL_UNDEF_R8

  do ll = 1, ng
     if (ng == 1) then
        if (EASE) then
           KK = iTable(:,4)
        else
           KK =[(k, k=1,ip)]
        endif
     else
        KK = iTable(:,5+ll)
     endif

     II = iTable(:,ll*2    )
     JJ = iTable(:,ll*2 + 1)

     where( rTable(:,3+ll) /=0.0)
        fr = rTable(:,3)/rTable(:,3+ll)
     endwhere

     if (ng == 1) then
       str_num=''
     else
       write(str_num, '(i0)') ll
     endif

     call formatter%put_var('i_indg'    //trim(str_num), II, rc=status)
     call formatter%put_var('j_indg'    //trim(str_num), JJ, rc=status)
     call formatter%put_var('frac_cell' //trim(str_num), fr, rc=status)
     call formatter%put_var('pfaf_index'//trim(str_num), KK, rc=status)
  enddo

  call formatter%put_var('min_lon', rTable(:, 6), rc=status)
  call formatter%put_var('max_lon', rTable(:, 7), rc=status)
  call formatter%put_var('min_lat', rTable(:, 8), rc=status)
  call formatter%put_var('max_lat', rTable(:, 9), rc=status)
  call formatter%put_var('elev',    rTable(:,10), rc=status)

  call formatter%close(rc=status)

  if (present(rc)) rc = status

end subroutine WriteTilingNC4

subroutine ReadTilingNC4(File, GridName, im, jm, nx, ny, iTable, rTable, srtm, rc)
  character(*),           intent(IN) :: File
  character(*), optional, intent(out) :: GridName(:)
  integer,      optional, intent(out) :: IM(:), JM(:)
  integer,      optional, intent(out) :: nx,ny
  integer,      optional, allocatable, intent(out) :: iTable(:,:)
  real(kind=8), optional, allocatable, intent(out) :: rTable(:,:)
  integer,      optional, intent(out) :: srtm
  integer,      optional, intent(out) :: rc

  type (Attribute), pointer     :: ref
  character(len=:), allocatable :: attr
  type (NetCDF4_FileFormatter)  :: formatter
  type (FileMetadata)           :: meta
  character(len=1)              :: str_num
  integer                       :: ng, ntile, status, ll
  class(*), pointer             :: attr_val(:)
  class(*), pointer             :: char_val
  integer, allocatable          :: tmp_int(:)
  real(kind=8), allocatable     :: fr(:)

  call formatter%open(File, pFIO_READ, rc=status)
  meta = formatter%read(rc=status)

  ng = 1
  if (meta%has_attribute("Grid1_Name")) ng = 2
  ntile = meta%get_dimension('tile')

  if (present(nx)) then
     ref => meta%get_attribute('raster_nx')
     attr_val => ref%get_values()
     select type(attr_val)
     type is (integer(INT32))
        nx = attr_val(1)
     endselect
  endif
  if (present(ny)) then 
     ref => meta%get_attribute('raster_ny')
     attr_val => ref%get_values()
     select type (attr_val)
     type is (integer(INT32))
        ny = attr_val(1)
     endselect
  endif

  if (present(srtm)) then
     ref => meta%get_attribute('SRTM_maxcat')
     attr_val => ref%get_values()
     select type (attr_val)
     type is (integer(INT32))
        srtm = attr_val(1)
     endselect
  endif

  do ll = 1, ng
    if (ng == 1) then
      str_num = ''
    else
      write(str_num, '(i0)') ll
    endif

    if (present(GridName)) then
       attr = 'Grid'//trim(str_num)//'_Name'
       ref =>meta%get_attribute(attr)
       char_val => ref%get_value()
       select type(char_val)
       type is(character(*))
          GridName(ll) = char_val
       class default
          print('unsupported subclass (not string) of attribute named '//attr)
       end select
    endif
    if (present(IM)) then
       attr = 'IM'//trim(str_num)
       ref =>meta%get_attribute(attr)
       attr_val => ref%get_values()
       select type(attr_val)
       type is( integer(INT32))
          IM(ll) = attr_val(1)
       end select
    endif
    if (present(JM)) then
       attr = 'JM'//trim(str_num)
       ref =>meta%get_attribute(attr)
       attr_val => ref%get_values()
       select type(attr_val)
       type is(integer(INT32))
          JM(ll) = attr_val(1)
       end select
    endif
  enddo

  if (present(iTable)) then
     allocate(iTable(ntile,0:7))
     allocate(tmp_int(ntile))
     call formatter%get_var('typ', iTable(:,0))
     do ll = 1, ng
       if (ng == 1) then
         str_num=''
       else
         write(str_num, '(i0)') ll
       endif

       call formatter%get_var('i_indg'    //trim(str_num), tmp_int, rc=status)
       iTable(:,ll*2) = tmp_int
       call formatter%get_var('j_indg'    //trim(str_num), tmp_int, rc=status)
       iTable(:,ll*2+1) = tmp_int
       call formatter%get_var('pfaf_index'//trim(str_num), tmp_int, rc=status)
       if ( ng == 1) then
         iTable(:,4) = tmp_int
       else
         iTable(:,5+ll) = tmp_int
       endif
     enddo
  endif

  if (present(rTable)) then
    allocate(rTable(ntile,10))
    allocate(fr(ntile))
    call formatter%get_var('com_lon', rTable(:,1),   rc=status)
    call formatter%get_var('com_lat', rTable(:,2),   rc=status)
    call formatter%get_var('area',    rTable(:,3),   rc=status)
    do ll = 1, ng
      if (ng == 1) then
        str_num=''
      else
        write(str_num, '(i0)') ll
      endif
      call formatter%get_var('frac_cell' //trim(str_num), fr, rc=status)
      where (fr /=0.0) rTable(:,3+ll) = rTable(:,3)/fr
    enddo
    call formatter%get_var('min_lon', rTable(:, 6), rc=status)
    call formatter%get_var('max_lon', rTable(:, 7), rc=status)
    call formatter%get_var('min_lat', rTable(:, 8), rc=status)
    call formatter%get_var('max_lat', rTable(:, 9), rc=status)
    call formatter%get_var('elev',    rTable(:,10), rc=status)
  endif

  if (present(rc)) rc= status
end subroutine

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

subroutine WriteLine(File, Unit, iTable, rTable, k, Zip, Verb)
  character*(*),     intent(IN) :: File
  integer,           intent(IN) :: Unit, k
  integer,           intent(IN) :: iTable(0:)
  real(kind=8),             intent(IN) :: rTable(:)
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
  real(kind=8)   :: fr
  real(kind=8)   :: xc, yc, area

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
  real(kind=8)   :: sphere, error
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


