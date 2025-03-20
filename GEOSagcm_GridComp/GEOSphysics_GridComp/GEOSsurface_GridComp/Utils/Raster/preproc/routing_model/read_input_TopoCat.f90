program main

  implicit none

  integer,parameter :: no=1459201
  integer,parameter :: nvl=3409
  integer,parameter :: nvo=3917
  integer,parameter :: nl=1426967

  integer,allocatable,dimension(:) :: lakeid_out,outid_out,lakeid_lake,lakeid_outV,outid_outV
  real,allocatable,dimension(:) :: lat_out,lon_out,lat_outV,lon_outV,&
                                 lakeaca_lake,lakearea_lake,lakeaca_outV,lakearea_outV

  real,allocatable,dimension(:) :: area_lake, aca_lake
  integer,allocatable,dimension(:) :: id_lake
  character(len=256) :: area_lake_file, id_lake_file, aca_lake_file
  integer :: i,j,k

  allocate(area_lake(nl),aca_lake(nl),id_lake(nl))
  ! Initialize file names
  area_lake_file = "input/Lake_area.csv"
  id_lake_file = "input/Hylak_id_lake.csv"
  aca_lake_file = "input/Cat_a_lake.csv"

  ! Read input data (You can implement your own read procedure if needed)
  open(77, file=area_lake_file, status="old")
  read(77, *) area_lake
  open(77, file=id_lake_file, status="old")
  read(77, *) id_lake
  open(77, file=aca_lake_file, status="old")
  read(77, *) aca_lake  

  ! Allocate arrays for filtered data
  allocate(lakearea_lake(nvl))
  allocate(lakeid_lake(nvl))
  allocate(lakeaca_lake(nvl))

  k = 0
  ! Filter lakes with area >= 50
  do i = 1, nl
     if (area_lake(i) .ge. 50.0) then
        k = k + 1
        lakearea_lake(k) = area_lake(i)
        lakeid_lake(k) = id_lake(i)
        lakeaca_lake(k) = aca_lake(i)
     end if
  end do
!-------------------------------------------------------------------------------------

  allocate(lakeid_out(no),outid_out(no),lat_out(no),lon_out(no))
  allocate(lakeid_outV(nvo),outid_outV(nvo),lat_outV(nvo),lon_outV(nvo),lakeaca_outV(nvo),lakearea_outV(nvo))

  open(77,file="input/Hylak_id_outlet.csv")
  read(77,*)lakeid_out
  open(77,file="input/Outlet_id.csv")
  read(77,*)outid_out
  open(77,file="input/Outlet_lat.csv")
  read(77,*)lat_out
  open(77,file="input/Outlet_lon.csv")
  read(77,*)lon_out

  k=0
  do i=1,no
    do j=1,nvl
  	  if(lakeid_out(i)==lakeid_lake(j))then
  	    k=k+1	
        outid_outV(k)=outid_out(i)
        lat_outV(k)=lat_out(i)
        lon_outV(k)=lon_out(i)
        lakeid_outV(k)=lakeid_lake(j)
        lakeaca_outV(k)=lakeaca_lake(j)
        lakearea_outV(k)=lakearea_lake(j)
  	  endif
    enddo
  enddo

  !print *,k

  open(88,file="temp/outlet_lat.txt")!
  do i=1,nvo
    write(88,*)lat_outV(i)
  enddo
  open(88,file="temp/outlet_lon.txt")!
  do i=1,nvo
    write(88,*)lon_outV(i)
  enddo
  open(88,file="temp/outlet_lakeid.txt")!
  do i=1,nvo
    write(88,*)lakeid_outV(i)
  enddo
  open(88,file="temp/outlet_lakeacaOBS.txt")!
  do i=1,nvo
    write(88,*)lakeaca_outV(i)
  enddo
  open(88,file="output/lake_outlet_lakearea.txt")!
  do i=1,nvo
    write(88,*)lakearea_outV(i)
  enddo


end