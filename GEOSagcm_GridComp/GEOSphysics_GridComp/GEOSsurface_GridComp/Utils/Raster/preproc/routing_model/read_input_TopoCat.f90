program main
!Main purpose: Reads lake and lake outlets information from Lake-TopoCat database.

use constant, only : no, nvl, nvo, nl=>nl_lake

  implicit none

  ! Declare arrays for outlet and lake data:
  integer, allocatable, dimension(:) :: lakeid_out, outid_out, lakeid_lake, lakeid_outV, outid_outV
  real, allocatable, dimension(:)    :: lat_out, lon_out, lat_outV, lon_outV, & 
                                        lakeaca_lake, lakearea_lake, lakeaca_outV, lakearea_outV

  ! Arrays for raw lake data:
  real, allocatable, dimension(:)    :: area_lake, aca_lake
  integer, allocatable, dimension(:) :: id_lake
  character(len=900)                 :: area_lake_file, id_lake_file, aca_lake_file
  integer                            :: i, j, k

  character(len=900)                 :: file_lake_area
  character(len=900)                 :: file_lake_id
  character(len=900)                 :: file_lake_aca
  character(len=900)                 :: file_lakeo_lakeid
  character(len=900)                 :: file_lakeo_id
  character(len=900)                 :: file_lakeo_lat
  character(len=900)                 :: file_lakeo_lon

  if (command_argument_count() /= 7) then
      print *, "no appropriate files found"
      stop
  endif
  call get_command_argument(1, file_lake_area)
  call get_command_argument(2, file_lake_id)
  call get_command_argument(3, file_lake_aca)
  call get_command_argument(4, file_lakeo_lakeid)
  call get_command_argument(5, file_lakeo_id)
  call get_command_argument(6, file_lakeo_lat)
  call get_command_argument(7, file_lakeo_lon)

  ! Allocate arrays for raw lake data (size nl):
  allocate(area_lake(nl), aca_lake(nl), id_lake(nl))
  
  ! Initialize file names for input lake data:
  area_lake_file = file_lake_area
  id_lake_file = file_lake_id
  aca_lake_file = file_lake_aca

  ! Read lake area, lake ID, and lake "aca" data from the input CSV files:
  open(77, file=trim(area_lake_file), status="old")
  read(77, *) area_lake
  open(77, file=trim(id_lake_file), status="old")
  read(77, *) id_lake
  open(77, file=trim(aca_lake_file), status="old")
  read(77, *) aca_lake  

  ! Allocate arrays for filtered lake data (size nvl):
  allocate(lakearea_lake(nvl))
  allocate(lakeid_lake(nvl))
  allocate(lakeaca_lake(nvl))

  k = 0
  ! Filter lakes: select only those with an area greater than or equal to 50.
  do i = 1, nl
     if (area_lake(i) .ge. 50.0) then
        k = k + 1
        lakearea_lake(k) = area_lake(i)
        lakeid_lake(k) = id_lake(i)
        lakeaca_lake(k) = aca_lake(i)
     end if
  end do
!-------------------------------------------------------------------------------------

  ! Allocate arrays for outlet data (raw arrays with size 'no'):
  allocate(lakeid_out(no), outid_out(no), lat_out(no), lon_out(no))
  ! Allocate arrays for matched outlet data (size 'nvo'):
  allocate(lakeid_outV(nvo), outid_outV(nvo), lat_outV(nvo), lon_outV(nvo), lakeaca_outV(nvo), lakearea_outV(nvo))

  ! Read outlet data from CSV files:
  open(77, file=trim(file_lakeo_lakeid))
  read(77, *) lakeid_out
  open(77, file=trim(file_lakeo_id))
  read(77, *) outid_out
  open(77, file=trim(file_lakeo_lat))
  read(77, *) lat_out
  open(77, file=trim(file_lakeo_lon))
  read(77, *) lon_out

  ! Match outlet records to filtered lakes:
  k = 0
  do i = 1, no
    do j = 1, nvl
      if (lakeid_out(i) == lakeid_lake(j)) then
        k = k + 1
        outid_outV(k) = outid_out(i)
        lat_outV(k) = lat_out(i)
        lon_outV(k) = lon_out(i)
        lakeid_outV(k) = lakeid_lake(j)
        lakeaca_outV(k) = lakeaca_lake(j)
        lakearea_outV(k) = lakearea_lake(j)
      end if
    end do
  end do

  ! Write matched outlet data to output text files:
  open(88, file="temp/outlet_lat.txt")
  do i = 1, nvo
    write(88, *) lat_outV(i)
  end do

  open(88, file="temp/outlet_lon.txt")
  do i = 1, nvo
    write(88, *) lon_outV(i)
  end do

  open(88, file="temp/outlet_lakeid.txt")
  do i = 1, nvo
    write(88, *) lakeid_outV(i)
  end do

  open(88, file="temp/outlet_lakeacaOBS.txt")
  do i = 1, nvo
    write(88, *) lakeaca_outV(i)
  end do

  open(88, file="output/lake_outlet_lakearea.txt")
  do i = 1, nvo
    write(88, *) lakearea_outV(i)
  end do

end