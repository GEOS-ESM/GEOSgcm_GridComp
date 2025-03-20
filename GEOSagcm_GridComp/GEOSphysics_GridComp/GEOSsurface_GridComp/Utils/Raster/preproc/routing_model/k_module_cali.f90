module k_module

  use river_read
   
  implicit none
  private
  public :: read_usgs_data, process_usgs_data, find_nearest_coords, get_station_inf, regression
  public :: filter_station, cal_Kmodel, get_valide_stations_gageii

contains
!------------------------------------------------------------
  subroutine read_usgs_data(nl, data)
    ! Subroutine to read USGS velocity and discharge data and store it in a 2D array
    integer, intent(in) :: nl
    real, allocatable, intent(out) :: data(:,:)
    
    character(len=100) :: var(2)
    character(len=256) :: filename
    character(len=100) :: line
    character(len=100) :: x(100)

    integer :: i, j, l, io, k
    integer, allocatable :: nv(:)

    !---------- Define the data types to be read ----------
    var = (/ "velocity", "discharge" /)

    ! Allocate the data array
    allocate(data(nl, 2))

    ! Loop over velocity and discharge files
    do l = 1, 2
      filename = "input/" // trim(var(l)) // ".txt"
      open(unit=77, file=trim(filename), status='old')

      ! Allocate temporary array for counting valid numbers per line
      allocate(nv(nl))
      
      ! Read each line from the file
      do k = 1, nl
        read(77, '(A)', iostat=io) line
        if (io /= 0) then
          print *, "Error reading line ", k, " from file: ", trim(filename)
          exit
        endif
        
        ! Read tokens from the line and store in array x
        do i = 1, 100
          read(line, *, iostat=io) (x(j), j=1, i)
          if (io == -1) then
            exit
          endif
        end do
        
        ! Count the number of valid values in the line
        nv(k) = i - 1
        
        ! If valid data exists, read the first value into the data array
        if (nv(k) >= 1) then
          read(x(1), *, iostat=io) data(k, l)
        else
          data(k, l) = -9999  ! Assign missing value if no data is available
        end if
      end do

      ! Deallocate the temporary array for valid number counts
      deallocate(nv)

      ! Close the file
      close(77)
    end do
  end subroutine read_usgs_data

!------------------------------------------------------------
  subroutine process_usgs_data(nl, ns, data, nv, nt, vel, dis)
    integer, intent(in) :: nl
    integer, intent(out) :: ns
    real, intent(inout),allocatable :: data(:,:)
    real, allocatable, intent(out) :: vel(:), dis(:)
    integer, allocatable,intent(out) :: nv(:)
    integer,intent(out) :: nt

    character(len=20), allocatable :: id(:)
    integer, allocatable :: nu(:)
    character(len=20), allocatable :: idu(:)
    integer :: i, k, ii

    ! Allocate arrays
    allocate(id(nl))


    ! Read IDs from the file
    open(unit=11, file="input/USGSID.txt", status="old")
    read(11, *) id
    close(11)

    ! Convert velocity from ft/s to m/s and discharge from ft^3/s to m^3/s
    !data(:,1) = data(:,1) * 0.3048     ! Convert ft/s to m/s
    !data(:,2) = data(:,2) * 0.0283168  ! Convert ft^3/s to m^3/s

    ! Initialize arrays
    k = 1
    ! Process ID and count occurrences of unique IDs
    do i = 2, nl
      if (.not.(trim(id(i)) == trim(id(i-1)))) then
        k = k + 1
      end if
    end do
    !print *, 'Number of unique IDs:', k
    ns=k
    allocate(nu(ns),nv(ns))
    allocate(idu(ns))
    ! Initialize arrays
    nu(1) = 1
    idu(1) = id(1)
    k = 1
    ! Process ID and count occurrences of unique IDs
    do i = 2, nl
      if (trim(id(i)) == trim(id(i-1))) then
        nu(k) = nu(k) + 1
      else
        k = k + 1
        nu(k) = 1
        idu(k) = id(i)
      end if
    end do

    ! Write idu to file (IDs without commas)
    open(unit=13, file="temp/id_for_site.txt")
    do i = 1, ns
      write(13, '(A)') trim(idu(i))//","
      !write(13, '(A)') trim(idu(i))
    end do
    close(13)
    open(unit=13, file="temp/id_for_site_nocomma.txt")
    do i = 1, ns
      write(13, '(A)') trim(idu(i))
      !write(13, '(A)') trim(idu(i))
    end do
    close(13)    

    ! Initialize variables
    nv = 0
    nv(1) = 1  
    k = 1
    ii = 0
    k = k - 1
    ! Assuming nm is defined and idA, vel, dis are allocated and filled
    ! Loop through the elements
    do i = 2, nl
      if (id(i) == id(i - 1)) then
        k = k + 1
        if (data(i,1) <= 0.0) then
          k = k - 1
        else if (data(i,2) <= 0.0) then
          k = k - 1
        end if
      else
        nv(ii + 1) = k  
        k = 1
        ii = ii + 1
        if (data(i,1) <= 0.0) then
          k = k - 1
        else if (data(i,2) <= 0.0) then
          k = k - 1
        end if
      end if
    end do
    nv(ii + 1) = k  
    !print *,"number valid records",sum(nv)

    nt=sum(nv)
    allocate(vel(nt),dis(nt))
    k=0
    do i=1,nl
      if(data(i,1)>0..and.data(i,2)>0.)then
        k=k+1
        vel(k)=data(i,1)
        dis(k)=data(i,2)
      endif
    enddo
    !open(unit=13, file="All_vel_trim.txt")
    !do i = 1, k
    !  write(13, *) vel(i)
    !end do  
    !open(unit=13, file="All_dis_trim.txt")
    !do i = 1, k
    !  write(13, *) dis(i)
    !end do


    ! Deallocate arrays
    deallocate(id)
    deallocate(nu)
    deallocate(idu)
    deallocate(data)

  end subroutine process_usgs_data  


!------------------------------------------------------------

  ! Subroutine to find the nearest latitude and longitude for each station
  subroutine find_nearest_coords(ns, nlat, nlon, lats, lons, lati, loni)
    integer, intent(in) :: ns, nlat, nlon
    real,allocatable,intent(inout) :: lats(:), lons(:)
    integer, allocatable, intent(out) :: lati(:), loni(:)

    real, allocatable :: lat1m(:), lon1m(:)
    real :: min_dist_lat, min_dist_lon, dist
    integer :: i, j, idx_min_lat, idx_min_lon

    ! Allocate arrays for lat/lon data
    allocate(lati(ns), loni(ns))
    allocate(lats(ns), lons(ns))
    allocate(lat1m(nlat), lon1m(nlon))

    !---- Read latitudes and longitudes for sites ----
    open(unit=10, file="input/lat_for_site_200.txt", status='old')
    do i = 1, ns
       read(10, *) lats(i)
    end do
    close(10)

    open(unit=11, file="input/lon_for_site_200.txt", status='old')
    do i = 1, ns
       read(11, *) lons(i)
    end do
    close(11)

    !---- Read 1-minute resolution lat/lon grid ----
    open(unit=12, file="input/lat_1m.txt", status='old')
    do i = 1, nlat
       read(12, *) lat1m(i)
    end do
    close(12)

    open(unit=13, file="input/lon_1m.txt", status='old')
    do i = 1, nlon
       read(13, *) lon1m(i)
    end do
    close(13)

    !---- Find nearest coordinates for each station ----
    do i = 1, ns
       ! Initialize minimum distance
       min_dist_lat = 1.0e20
       min_dist_lon = 1.0e20
       idx_min_lat = -1
       idx_min_lon = -1

       ! Find nearest latitude
       do j = 1, nlat
          dist = abs(lats(i) - lat1m(j))
          if (dist < min_dist_lat) then
             min_dist_lat = dist
             idx_min_lat = j
          end if
       end do
       lati(i) = idx_min_lat

       ! Find nearest longitude
       do j = 1, nlon
          dist = abs(lons(i) - lon1m(j))
          if (dist < min_dist_lon) then
             min_dist_lon = dist
             idx_min_lon = j
          end if
       end do
       loni(i) = idx_min_lon
    end do

   !---- Write output files with Fortran indexing (1-based) ----
    !open(unit=14, file="USGS_data/lati_for_site_200.txt", status='replace')
    !do i = 1, ns
    !   write(14, *) lati(i)
    !end do
    !close(14)

    !open(unit=15, file="USGS_data/loni_for_site_200.txt", status='replace')
    !do i = 1, ns
    !   write(15, *) loni(i)
    !end do
    !close(15)

    ! Deallocate arrays
    deallocate(lat1m)
    deallocate(lon1m)
  end subroutine find_nearest_coords
!------------------------------------------------------------
  subroutine get_station_inf(ns, nc, nlat, nlon, lati, loni, catid, Qclmt, slp, KImodel_all,exp_slp,exp_clmt,fac_str)
    integer, intent(in) :: ns, nc, nlat, nlon
    integer, intent(in) :: lati(nlon), loni(nlon)
    integer, allocatable, intent(out) :: catid(:)
    real,allocatable,intent(out) :: Qclmt(:),slp(:)
    real,allocatable,intent(out) :: KImodel_all(:)
    real,intent(in) :: exp_slp,exp_clmt,fac_str

    integer,allocatable :: catchind(:,:)
    real,allocatable,dimension(:) :: Qclmt_all,slp_all,Kstr_all,Qstr_all
    integer :: i

    allocate(catchind(nlon,nlat),catid(ns))
    allocate(Qclmt_all(nc),slp_all(nc))
    allocate(Qclmt(ns),slp(ns))
    allocate(KImodel_all(nc),Kstr_all(nc),Qstr_all(nc))    
    
    call read_ncfile_int2d("input/SRTM_PfafData.nc","CatchIndex",catchind,nlon,nlat)

    do i=1,ns
      catid(i)=catchind(loni(i),lati(i))
    end do

    open(88,file="temp/catid_for_site_200.txt")
    do i=1,ns
      write(88,*)catid(i)
    end do    

    open(77,file="output/Pfaf_qri.txt")
    read(77,*)Qclmt_all
    where(Qclmt_all<1.e-8) Qclmt_all=1.e-8
    open(77,file="temp/Pfaf_slope.txt")
    read(77,*)slp_all           
    open(77,file="output/Pfaf_qstr.txt")
    read(77,*)Qstr_all    
    where(Qstr_all<1.e-8) Qstr_all=1.e-8

    do i=1,ns
      if(catid(i)/=-9999)then
        Qclmt(i)=Qclmt_all(catid(i))
        slp(i)=slp_all(catid(i))
      else
        Qclmt(i)=-9999
        slp(i)=-9999
      endif
    enddo 

    KImodel_all = (Qclmt_all**(exp_clmt)) * (slp_all**(exp_slp))

    Kstr_all = fac_str * (Qstr_all**(exp_clmt)) * (slp_all**(exp_slp))

    open(88,file="output/Pfaf_Kstr_PR_fac1_0p35_0p45_0p2_n0p2.txt")
    do i=1,nc
      write(88,*)Kstr_all(i)
    enddo    


    !open(88,file="USGS_data/qri_for_site_200.txt")
    !do i=1,ns
    !  write(88,*)Qclmt(i)
    !end do   
    !open(88,file="USGS_data/slp_for_site_200.txt")
    !do i=1,ns
    !  write(88,*)slp(i)
    !end do       
    
    
    deallocate(catchind,Qclmt_all,slp_all,Kstr_all,Qstr_all)

  end subroutine get_station_inf
!------------------------------------------------------------
  subroutine get_valide_stations_gageii(ns,nc,catid_sta,flag_thres)
    integer,intent(in) :: ns,nc
    integer,intent(in) :: catid_sta(ns)
    integer,allocatable,intent(out) :: flag_thres(:)


  integer, parameter :: nga = 9067
  integer, parameter :: nv = 5704  
  real,    parameter :: thr_sel = 0.3

  real,dimension(:),allocatable :: acar_pfaf

  integer :: i, j, k, cid
  character(len=20) :: id_gages(nga)
  character(len=20) :: id_sta(ns)
  integer           :: flag_gageii(ns)
  real :: acar_gages(nga)
  real :: acar_gages_sta(ns),acar_sta(ns)
  character(len=20) :: line
  integer :: ios

  allocate(flag_thres(ns))
  
  ! Initialize acar_6156 array with missing value
  acar_sta = -9999.0
  k = 0
  
  ! Read id_gages from file
  open(unit=10, file="input/id_gagesii.txt", status="old", action="read")
  do j = 1, nga
    read(10,'(A)', iostat=ios) id_gages(j)
    if (ios /= 0) then
      print *, "Error reading id_gagesii.txt"
      stop
    end if
  end do
  close(10)
  
  ! Read id_6156 from file
  open(unit=11, file="temp/id_for_site_nocomma.txt", status="old", action="read")
  do i = 1, ns
    read(11,'(A)', iostat=ios) id_sta(i)
    if (ios /= 0) then
      print *, "Error reading id_for_site_nocomma.txt"
      stop
    end if
  end do
  close(11)
  
  ! Read acar_gages from file
  open(unit=12, file="input/acar_gagesii.txt", status="old", action="read")
  do j = 1, nga
    read(12,*, iostat=ios) acar_gages(j)
    if (ios /= 0) then
      print *, "Error reading acar_gagesii.txt"
      stop
    end if
  end do
  close(12)
  
  flag_gageii = 0
  ! Compare id_sta and id_gages, and update acar_sta if there's a match
  do i = 1, ns
    do j = 1, nga
      if (trim(id_gages(j)) == trim(id_sta(i))) then
        acar_gages_sta(i) = acar_gages(j)
        flag_gageii(i) = 1
        k = k + 1
        exit  ! Exit inner loop if match is found
      end if
    end do
  end do

  print *, "Number of matches:", sum(flag_gageii)

  allocate(acar_pfaf(nc))
  open(77,file="temp/Pfaf_acar.txt")
  read(77,*)acar_pfaf

  do i = 1, ns
    if(catid_sta(i)/=-9999)then
      cid = catid_sta(i)
      acar_sta(i) = acar_pfaf(cid)
    else
      acar_sta(i) = -9999.
    end if
  end do  


  flag_thres = 0
  do i = 1, ns
    if(flag_gageii(i)==1 .and. catid_sta(i)/=-9999)then
      if(acar_sta(i).ge.(1.-thr_sel)*acar_gages_sta(i) .and. acar_sta(i).le.(1.+thr_sel)*acar_gages_sta(i))then
        flag_thres(i) = 1
      endif
    endif
  end do

  print *,"Number of valid:", sum(flag_thres)

  deallocate(acar_pfaf)
  !open(88,file="flag_thr03_7065.txt")
  !do i = 1,ns
  !  write(88,*)flag_thres(i)
  !enddo

  end subroutine get_valide_stations_gageii
!------------------------------------------------------------
  subroutine regression(nt,vel_ori,dis_ori,nv,ns,Qclmt,slp,KKobs,KImodel,exp_slp,exp_clmt,mm,MU)
    integer,intent(in) :: nt, ns
    real,intent(inout),allocatable :: vel_ori(:), dis_ori(:)
    integer,intent(in) :: nv(ns)
    real,intent(inout),allocatable :: Qclmt(:),slp(:)
    real,intent(out),allocatable :: KKobs(:),KImodel(:)
    real,intent(in) :: exp_slp,exp_clmt,mm,MU

    real,allocatable,dimension(:) :: x,y,yest
 
    integer :: thres=100
    integer :: i,j
    real :: k(ns),cdtm(ns),med
    integer :: acc(ns)
    real,allocatable :: vel(:), dis(:)

    allocate(vel(nt),dis(nt))
    vel=vel_ori*0.3048 !m/s
    dis=dis_ori*0.0283168 !m3/s

    acc(1)=nv(1)
    do i=2,ns
      acc(i)=acc(i-1)+nv(i)
    end do
    !open(88,file="USGS_data/acc_noMISSING_200.txt")
    !do i=1,ns
    !  write(88,*)acc(i)
    !end do    
    !print *,"5.1"
    do i=1,ns
      if(nv(i)>=thres)then
        allocate( x(nv(i)), y(nv(i)), yest(nv(i))) 
        x=dis( acc(i)-nv(i)+1:acc(i) )**mm
        y=vel( acc(i)-nv(i)+1:acc(i) )
        k(i)=sum(x*y)/sum(x*x)
        yest=k(i)*x
        cdtm(i)=cal_cdtm(y,yest)
        deallocate(x,y,yest)
      else
        k(i)=-9999.
        cdtm(i)=-9999.
      endif
    enddo
    med=median(cdtm)

    where(cdtm<0.5)k=-9999.
    !print *,"mm=",mm,",cdtm_med=",med,",stop now!"

    !print *,"5.2"
    allocate(KKobs(ns))
    do i=1,ns
      if(k(i)/=-9999.and.Qclmt(i)/=-9999.)then
        KKobs(i)=k(i)/(Qclmt(i)**(MU-mm))
      else
        KKobs(i)=-9999.
      endif
    end do
    
    !open(88,file="KKobs_mm0p40_MU0p10_7065.txt")
    !do i=1,ns
    !  write(88,*)KKobs(i)
    !enddo
    
    !print *,"mm=",mm,",cdtm_med=",med,",stop now!"
    !stop    

    !open(88,file="USGS_data/KKobs_200.txt")
    !do i=1,ns
    !  write(88,*)KKobs(i)
    !end do  

    allocate(KImodel(ns))
    KImodel = (Qclmt**(exp_clmt)) * (slp**(exp_slp))

    deallocate(vel,dis)
    !deallocate(vel_ori,dis_ori) 

  end subroutine regression
!------------------------------------------------------------
  subroutine filter_station(nc,ns,np,lats_full,lons_full,Qclmt_full,slp_full,catid_full,KKobs_full,KImodel_full,Qclmt,slp,catid,KKobs,KImodel,flag_gageii)
    integer,intent(in) :: ns,nc
    integer,intent(out) :: np
    real,intent(inout),allocatable :: lats_full(:),lons_full(:),Qclmt_full(:),slp_full(:),KKobs_full(:),KImodel_full(:)
    real,intent(out),allocatable :: Qclmt(:),slp(:),KKobs(:),KImodel(:)
    integer,intent(inout),allocatable :: catid_full(:)
    integer,intent(out),allocatable :: catid(:)  
    integer,intent(inout),allocatable :: flag_gageii(:)

    integer,allocatable :: flag_slp(:)
    real,allocatable :: lats(:),lons(:)
    integer :: i,k
    integer,allocatable :: flag_7065(:)


    allocate(flag_slp(nc))
    open(77,file="temp/Pfaf_slope_flag.txt")
    read(77,*)flag_slp 

    allocate(flag_7065(ns))
    flag_7065=0

    !open(77,file="flag_thr03_7065.txt")
    !read(77,*)flag_gageii

    k=0
    do i=1,ns
      if(catid_full(i).ne.-9999.and.KKobs_full(i)/=-9999..and.slp_full(i)>1.e-5.and.flag_slp(catid_full(i))==1.and.flag_gageii(i)==1)then
!      if(catid_full(i).ne.-9999.and.KKobs_full(i)/=-9999..and.slp_full(i)>1.e-5.and.flag_gageii(i)==1)then
        k=k+1
      endif
    enddo
    np=k
    print *,"number of valid stations: ",np
    !stop
    allocate(Qclmt(np),slp(np),catid(np),KKobs(np),KImodel(np))
    allocate(lats(np),lons(np))
    k=0
    do i=1,ns
      if(catid_full(i).ne.-9999.and.KKobs_full(i)/=-9999..and.slp_full(i)>1.e-5.and.flag_slp(catid_full(i))==1.and.flag_gageii(i)==1)then
!      if(catid_full(i).ne.-9999.and.KKobs_full(i)/=-9999..and.slp_full(i)>1.e-5.and.flag_gageii(i)==1)then
        k=k+1
        Qclmt(k)=Qclmt_full(i)
        slp(k)=slp_full(i)
        KKobs(k)=KKobs_full(i)
        KImodel(k)=KImodel_full(i)
        catid(k)=catid_full(i)
        lats(k)=lats_full(i)
        lons(k)=lons_full(i)
        flag_7065(i)=1
      endif
    enddo    
    !open(88,file="flag_7065_stations_1265.txt")
    !do i=1,ns
    !  write(88,*)flag_7065(i)
    !enddo
    !stop    
    !open(88,file="lats_stations.txt")
    !do i=1,np
    !  write(88,*)lats(i)
    !enddo
    !open(88,file="lons_stations.txt")
    !do i=1,np
    !  write(88,*)lons(i)
    !enddo    

    deallocate(Qclmt_full,slp_full,KKobs_full,KImodel_full,flag_slp,flag_gageii,lats,lons)

  end subroutine filter_station
!------------------------------------------------------------
  subroutine cal_Kmodel(ns,np,nc,MU,exp_slp,exp_clmt,Qclmt,slp,KKobs,KImodel,KImodel_all,catid,catid_full,ccr,rms)
    integer,intent(in) :: ns,np,nc
    real,intent(in) :: MU,exp_slp,exp_clmt
    real,intent(inout),allocatable :: Qclmt(:),slp(:),KKobs(:),KImodel(:)
    real,intent(inout),allocatable :: KImodel_all(:)
    integer,intent(inout),allocatable :: catid(:),catid_full(:)
    real,intent(inout) :: ccr,rms
    
    real,allocatable :: KKobs_sort(:), KImodel_sort(:), KKmodel_full(:)
    real, allocatable, dimension(:) :: dis,sca,Kv,KKmodel
    integer,allocatable,dimension(:) :: gear

    character(len=50) :: MU_s,exp_slp_s,exp_clmt_s

    integer :: bulk,i,lev
    real :: Kper(11),KMper(11),rat(11),dis_full(11)

    write(MU_s,'(f4.2)')MU
    write(exp_slp_s,'(f4.2)')exp_slp
    if(exp_clmt>=0.)then
      write(exp_clmt_s,'(f4.2)')exp_clmt
    else
      write(exp_clmt_s,'(f4.2)') -1.*exp_clmt
      exp_clmt_s="n"//trim(exp_clmt_s)
    endif

    allocate(KKobs_sort(np),KImodel_sort(np))
    call sort(np,KKobs,KKobs_sort)
    call sort(np,KImodel,KImodel_sort)

    bulk=np/10
    Kper(1)=KKobs_sort(1)
    KMper(1)=KImodel_sort(1)
    do i=2,10
      Kper(i)=KKobs_sort(bulk*(i-1))
      KMper(i)=KImodel_sort(bulk*(i-1))
    enddo
    Kper(11)=KKobs_sort(np)
    KMper(11)=KImodel_sort(np)
    rat=Kper/KMper

    !open(88,file="rat_Kper2KMper.txt")
    !do i=1,11
    !  write(88,*)rat(i)
    !enddo
    !close(88)
    !exit

    allocate(gear(nc),dis(nc),sca(nc),Kv(nc))
    
    gear=12
    dis=-9999.
    do i=1,nc
      do lev=1,11
        if(KImodel_all(i)<=KMper(lev))then
          gear(i)=lev
          dis(i)=KMper(lev)-KImodel_all(i)
          exit
        endif
      end do
    enddo    

    dis_full(1)=KMper(1)
    do i=2,11
      dis_full(i)=KMper(i)-KMper(i-1)
    enddo

    do i=1,nc
      if(gear(i)==1)then
        sca(i)=rat(1)
      elseif(gear(i)==12)then
        sca(i)=rat(11)
      else
        sca(i)= ( rat(gear(i)-1)*dis(i) + rat(gear(i))*(dis_full(gear(i))-dis(i)) ) / dis_full(gear(i))
      endif 
      Kv(i)=KImodel_all(i)*sca(i)
    enddo    

    open(88,file="output/Pfaf_Kv_PR_0p35_0p45_0p2_n0p2.txt")
!    open(88,file="Pfaf_Kv_PR_0p4_0p1_0p5_0p2.txt")
    do i=1,nc
      write(88,*)Kv(i)
    enddo

    allocate(KKmodel_full(ns))
    do i=1,ns
      if(catid_full(i)/=-9999)then
        KKmodel_full(i)=Kv(catid_full(i))
      else
        KKmodel_full(i)=-9999.
      endif
    enddo

    !open(88,file="KKmodel_7065/KKmodel_7065_"//trim(MU_s)//"_"//trim(exp_slp_s)//"_"//trim(exp_clmt_s)//".txt")
    !do i=1,ns
    !  write(88,*)KKmodel_full(i)
    !enddo   

    allocate(KKmodel(np))
    do i=1,np
      KKmodel(i)=Kv(catid(i))
    enddo

    ccr=cal_ccr(KKobs,KKmodel)
    rms=cal_rms(KKobs,KKmodel,np)


    !open(88,file="KKobs_stations.txt")
    !do i=1,np
    !  write(88,*)KKobs(i)
    !enddo       
    !print *,ccr


    deallocate(KKobs_sort,KImodel_sort)
    deallocate(KImodel_all,gear,dis,sca,Kv)
    deallocate(Qclmt,slp,KKobs,KImodel,catid,KKmodel,catid_full,KKmodel_full)
  end subroutine cal_Kmodel

  subroutine sort(np, data, data_sort)
    integer, intent(in) :: np           ! The size of the array
    real, intent(in) :: data(np)        ! Input array to be sorted
    real, intent(out) :: data_sort(np)  ! Output sorted array
    integer :: i, j
    real :: temp

    ! Copy input array to output array
    data_sort = data

    ! Perform a bubble sort (simple sorting algorithm)
    do i = 1, np-1
      do j = 1, np-i
        if (data_sort(j) > data_sort(j+1)) then
          ! Swap the elements
          temp = data_sort(j)
          data_sort(j) = data_sort(j+1)
          data_sort(j+1) = temp
        end if
      end do
    end do
  end subroutine sort

  function cal_ccr(y, yest) result(ccr)
    real, intent(in) :: y(:)
    real, intent(in) :: yest(:)
    real :: ccr
    real :: mean_y, mean_yest
    real :: sum_y, sum_yest
    real :: sum_num, sum_den_y, sum_den_yest
    integer :: n
    integer :: i

    n = size(y)
    if (n /= size(yest)) then
      print *, "Error: Arrays must have the same length"
      ccr = 0.0
      return
    endif

    ! Calculate means
    sum_y = sum(y)
    sum_yest = sum(yest)
    mean_y = sum_y / n
    mean_yest = sum_yest / n

    ! Calculate numerator and denominators for correlation coefficient
    sum_num = 0.0
    sum_den_y = 0.0
    sum_den_yest = 0.0
    do i = 1, n
      sum_num = sum_num + (y(i) - mean_y) * (yest(i) - mean_yest)
      sum_den_y = sum_den_y + (y(i) - mean_y) ** 2
      sum_den_yest = sum_den_yest + (yest(i) - mean_yest) ** 2
    end do

    ! Calculate correlation coefficient
    if (sum_den_y == 0.0 .or. sum_den_yest == 0.0) then
      print *, "Error: Zero variance in input arrays"
      ccr = 0.0
    else
      ccr = sum_num / sqrt(sum_den_y * sum_den_yest)
    end if

  end function cal_ccr  

 function cal_rms(k_obs,k_model, n) result(rms)
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: k_obs(n),k_model(n)
  real :: rms
  real :: sum_sq_diff
  integer :: i

  sum_sq_diff = 0.0

  do i = 1, n
    sum_sq_diff = sum_sq_diff + ((k_model(i) - k_obs(i)) / k_obs(i))**2
  end do

  rms = sqrt(sum_sq_diff / n)
end function cal_rms 

  function cal_cdtm(y, yest) result(dtmc)
    real, intent(in) :: y(:)
    real, intent(in) :: yest(:)
    real :: dtmc
    real :: ss_tot, ss_res
    real :: mean_y
    integer :: n, i

    n = size(y)
    if (n /= size(yest)) then
      print *, "Error: Arrays must have the same length"
      dtmc = 0.0
      return
    endif

    ! Calculate mean of y
    mean_y = sum(y) / n

    ! Calculate total sum of squares (SS_tot)
    ss_tot = sum((y - mean_y)**2)

    ! Calculate residual sum of squares (SS_res)
    ss_res = sum((y - yest)**2)

    ! Calculate coefficient of determination (R^2)
    if (ss_tot == 0.0) then
      print *, "Error: Zero total sum of squares"
      dtmc = 0.0
    else
      dtmc = 1.0 - (ss_res / ss_tot)
    endif

  end function cal_cdtm

function median(data) result(med)
    implicit none
    real, intent(in) :: data(:)  
    real :: med                  
    real :: sorted_data(size(data))  
    integer :: n_valid           
    integer :: i                 
    
    n_valid = 0
    do i = 1, size(data)
        if (data(i) /= -9999.0) then
            n_valid = n_valid + 1
            sorted_data(n_valid) = data(i)
        end if
    end do

    if (n_valid == 0) then
        med = -9999.0  
        return
    end if

    call sort2(sorted_data(1:n_valid))
    
    if (mod(n_valid, 2) == 0) then
        med = (sorted_data(n_valid/2) + sorted_data(n_valid/2 + 1)) / 2.0
    else
        med = sorted_data((n_valid + 1) / 2)
    end if

end function median

subroutine sort2(arr)
    implicit none
    real, intent(inout) :: arr(:)
    integer :: i, j
    real :: temp
    
    do i = 2, size(arr)
        temp = arr(i)
        j = i - 1
        do while (j >= 1 .and. arr(j) > temp)
            arr(j + 1) = arr(j)
            j = j - 1
        end do
        arr(j + 1) = temp
    end do
end subroutine sort2

!------------------------------------------------------------

end module k_module