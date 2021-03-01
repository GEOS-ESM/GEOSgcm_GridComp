MODULE MOD_GEOS_IO

  USE netcdf
  
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'


contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE NC_READ_INPUT_FILE -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE NC_READ_INPUT_FILE(fname,Npnts,Nl,Nhydro,lon,lat,p,ph,z,zh,T,qv,rh,tca,cca, &
            mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_ccrain,fl_ccsnow,Reff, &
            dtau_s,dtau_c,dem_s,dem_c,skt,landmask,emsfc_lw,mode,Nlon,Nlat)
  
  !USE netcdf

  !IMPLICIT NONE



  !INCLUDE 'netcdf.inc'
    
    !Arguments
    character(len=512),intent(in) :: fname ! File name
    integer,intent(in) :: Npnts,Nl,Nhydro
    real,dimension(Npnts),intent(out) :: lon,lat
    real,dimension(Npnts,Nl),target,intent(out) :: p,ph,z,zh,T,qv,rh,tca,cca, &
                  mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_ccrain,fl_ccsnow, &
                  dtau_s,dtau_c,dem_s,dem_c
    real,dimension(Npnts,Nl,Nhydro),intent(out) :: Reff
    real,dimension(Npnts),intent(out) :: skt,landmask
    real,intent(out) :: emsfc_lw
    integer,intent(out) :: mode,Nlon,Nlat

    !Local variables
    integer :: Npoints,Nlevels,i,j,k
    character(len=128) :: vname
    integer,parameter :: NMAX_DIM=5
    integer :: vrank,vdimid(NMAX_DIM)
    character(len=MAXNCNAM) :: dimname(NMAX_DIM)
    integer :: ncid,vid,ndims,nvars,ngatts,recdim,dimsize(NMAX_DIM)
    integer :: errst
    logical :: Llat,Llon,Lpoint
    integer :: Na,Nb,Nc,Nd,Ne
    real,dimension(Npnts) :: ll
    integer,dimension(:),allocatable :: plon,plat
    real,allocatable :: x1(:),x2(:,:),x3(:,:,:),x4(:,:,:,:),x5(:,:,:,:,:) ! Temporary arrays



    
    ! Open file
    errst = nf90_open(fname, nf90_nowrite, ncid)
    
    ! Get information about dimensions. Curtain mode or lat/lon mode?
    Llat  =.false.
    Llon  =.false.
    Lpoint=.false.
    errst = nf90_inquire(ncid, ndims, nvars, ngatts, recdim)



   ! Close file
    call ncclos(ncid,errst)
  END SUBROUTINE NC_READ_INPUT_FILE
end MODULE MOD_GEOS_IO
