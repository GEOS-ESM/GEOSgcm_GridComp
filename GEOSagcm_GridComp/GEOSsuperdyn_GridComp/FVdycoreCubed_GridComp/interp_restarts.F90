#define VERIFY_(A)   if((A)/=0) then; PRINT *, 'Interp_restarts.x', __LINE__; call MPI_Abort(A); endif

program interp_restarts
#define DEALLOCGLOB_(A) if(associated(A)) then;A=0;call MAPL_DeAllocNodeArray(A,rc=STATUS);if(STATUS==MAPL_NoShm) deallocate(A,stat=STATUS);NULLIFY(A);endif

!--------------------------------------------------------------------!
! purpose: driver for interpolation of GEOS FV and Moist restarts    !
!          to the cubed-sphere grid with optional vertical levels    !
!--------------------------------------------------------------------!
   use ESMF
   use mpp_mod,        only: mpp_error, FATAL, NOTE, mpp_root_pe, mpp_broadcast
   use fms_mod,        only: print_memory_usage, fms_init, fms_end, file_exist
   use fv_control_mod, only: fv_init1, fv_init2, fv_end
   use fv_arrays_mod,  only: fv_atmos_type, REAL4, REAL8, FVPRC
   use fv_mp_mod,      only: is_master, ng, mp_gather, tile
   use fv_regrid_c2c, only: get_geos_ic
   use fv_regridding_utils
   use fv_grid_utils_mod, only: ptop_min
   use init_hydro_mod, only: p_var
   use constants_mod,  only: pi, omega, grav, kappa, rdgas, rvgas, cp_air
   use fv_diagnostics_mod,only: prt_maxmin
! use fv_eta_mod,     only: set_eta
   use m_set_eta,     only: set_eta
   use memutils_mod, only: print_memuse_stats
   use MAPL_IOMod
   use MAPL_ShmemMod
   use MAPL_ConstantsMod
   use rs_scaleMod
   use pFIO_StringVectorMod
   use pFIO_StringIntegerMapMod

   implicit none

#include "mpif.h"

   type(fv_atmos_type), allocatable, save :: FV_Atm(:)
   logical, allocatable, save             :: grids_on_this_pe(:)

   real, parameter:: zvir = rvgas/rdgas - 1.

   character(ESMF_MAXSTR) :: fname1, fname2, str, astr
#ifndef __GFORTRAN__
   external :: getarg, iargc
   integer iargc
#endif
   real(FVPRC) :: dt

   real(ESMF_KIND_R8), allocatable :: r8_ak(:)
   real(ESMF_KIND_R8), allocatable :: r8_bk(:)
   real(ESMF_KIND_R8), allocatable :: r8_akbk(:)

   real(ESMF_KIND_R4), pointer :: r4_local(:,:,:)
   real(ESMF_KIND_R8), pointer :: r8_local(:,:,:), pt_local(:,:,:)
   real(ESMF_KIND_R4), pointer :: r4_local2D(:,:)

   integer i,j,k,n,iq, itmp, ihydro
   integer im,jm,km,nq
   real(ESMF_KIND_R8) :: ptop
   real(ESMF_KIND_R8) :: pint

   integer :: is,ie, js,je
   integer :: ks
   integer :: status
   integer :: header(6)
   integer :: IUNIT=15
   !integer :: OUNIT=17
   integer :: ounit

   integer :: nmoist
   logical :: isBinFV, isBinMoist
   integer :: ftype
   type(MAPL_NCIO) :: ncio,ncioOut
   integer :: nVars,imc,jmc,lonid,latid,levid,edgeid
   character(62) :: vname
   type(StringVector) :: moist_variables
   type(StringIntegerMap) :: tracer_names
   character(len=128) :: moist_order(9) = (/"Q   ","QLLS","QLCN","CLLS","CLCN","QILS","QICN","NCPL","NCPI"/)
   integer :: p_split, npx, npy, npz,ivar,lcnt_var,iq0
   integer :: n_args,n_files,nlevs,nedges,ifile,nlev,n_output,nfv_vars
   character(len=ESMF_MAXPATHLEN), allocatable :: extra_files(:),extra_output(:)
   type(fv_rst), pointer :: rst_files(:) => null()
   type(ArrDescr) :: ArrDes
   integer        :: info
   logical        :: amWriter
   integer :: isl,iel,jsl,jel,npes_x,npes_y,n_writers,n_readers,ndims,dimsize(3)
   type(ESMF_Grid) :: grid
   logical :: scale_rst

   real(FVPRC), allocatable :: q_wat(:,:,:,:)

! Start up FMS/MPP
   print_memory_usage = .true.
   call fms_init()
   call ESMF_Initialize(logKindFlag=ESMF_LOGKIND_NONE,mpiCommunicator=MPI_COMM_WORLD)
   p_split = 1
   call fv_init1(FV_Atm, dt, grids_on_this_pe, p_split)
   call print_memuse_stats('interp_restarts: fms_init')

   n_args = command_argument_count()
   n_files = 0
   n_output = 0
   n_writers=1
   n_readers=1
   ihydro = 1
   scale_rst = .true.
   do i=1,n_args
     call get_command_argument(i,str)
     select case(trim(str))
     case ('-im')
        call get_command_argument(i+1,astr)
        read(astr,*)npx
     case('-lm')
        call get_command_argument(i+1,astr)
        read(astr,*)npz
     case('-do_hydro')
        call get_command_argument(i+1,astr)
        read(astr,*)ihydro
     case('-input_files')
        do j=i+1,n_args
           call get_command_argument(j,astr)
           if ( index(astr,'-') .ne. 0) then
              exit
           end if
           n_files=n_files+1
        enddo
        allocate(extra_files(n_files))
        nq = 0
        do j=i+1,i+n_files
           nq=nq+1
           call get_command_argument(j,extra_files(nq))
        enddo
     case('-output_files')
        do j=i+1,n_args
           call get_command_argument(j,astr)
           if ( index(astr,'-') .ne. 0) then
              exit
           end if
           n_output=n_output+1
        enddo
        allocate(extra_output(n_output))
        nq = 0
        do j=i+1,i+n_output
           nq=nq+1
           call get_command_argument(j,extra_output(nq))
        enddo
     case('-nreader')
        call get_command_argument(i+1,astr)
        read(astr,*)n_readers
     case('-nwriter')
        call get_command_argument(i+1,astr)
        read(astr,*)n_writers
     case('-scalers')
        call get_command_argument(i+1,astr)
        if (trim(astr) == "T") then
           scale_rst=.true.
        else if (trim(astr) == "F") then
           scale_rst=.false.
        else
           write(*,*)'bad argument to scalers, will scale by default'
        end if
     end select
   end do


   npx = npx+1
   FV_Atm(1)%flagstruct%npx=npx
   npy = npx
   FV_Atm(1)%flagstruct%npy=npy

   FV_Atm(1)%flagstruct%npz=npz
   FV_Atm(1)%flagstruct%ntiles = 6

   FV_Atm(1)%flagstruct%hydrostatic = .true.
   if (ihydro == 0) FV_Atm(1)%flagstruct%hydrostatic = .false.
   FV_Atm(1)%flagstruct%Make_NH = .false.
   if (.not. FV_Atm(1)%flagstruct%hydrostatic) FV_Atm(1)%flagstruct%Make_NH = .true.

   if (n_files > 0) allocate(rst_files(n_files)) 

! Initialize SHMEM in MAPL
   call MAPL_GetNodeInfo (comm=MPI_COMM_WORLD, rc=status)
   call MAPL_InitializeShmem (rc=status)

                       write(fv_atm(1)%flagstruct%grid_file, "('c',i2.2,'_mosaic.nc')") npx-1
   if (npx-1 >=   100) write(fv_atm(1)%flagstruct%grid_file, "('c',i3.3,'_mosaic.nc')") npx-1
   if (npx-1 >=  1000) write(fv_atm(1)%flagstruct%grid_file, "('c',i4.4,'_mosaic.nc')") npx-1
   if (npx-1 >= 10000) write(fv_atm(1)%flagstruct%grid_file, "('c',i5.5,'_mosaic.nc')") npx-1
   dt = 1800
   call fv_init2(FV_Atm, dt, grids_on_this_pe, p_split)

   if (size(extra_files) > 0) then
      if (size(extra_files) /= size(extra_output)) call mpp_error(FATAL, 'the number of extra input and output file names must be same size')
   end if
   call print_memuse_stats('interp_restarts: fv_init')

! Determine Total Number of Tracers (MOIST, GOCART, PCHEM, ANA)
! -------------------------------------------------------------
   nmoist  = 0
   isBinFV     = .true.
   isBinMoist  = .true.

   if (is_master()) print*, 'HYDROSTATIC : ', FV_Atm(1)%flagstruct%hydrostatic  
   if (is_master()) print*, 'Make_NH     : ', FV_Atm(1)%flagstruct%Make_NH
   if (is_master()) print*, 'Tracers     : ', FV_Atm(1)%ncnst

! Need to get ak/bk
   if( file_exist("fvcore_internal_restart_in") ) then
      call MAPL_NCIOGetFileType("fvcore_internal_restart_in",ftype)
      if (ftype /= 0) then
         isBinFV = .true.
      else 
         isBinFV = .false.
      end if
      if (isBinFV) then 
         open(IUNIT,file='fvcore_internal_restart_in' ,access='sequential',form='unformatted',status='old')
! Headers
         read (IUNIT, IOSTAT=status) header
         read (IUNIT, IOSTAT=status) header(1:5)
         im=header(1)
         jm=header(2)
         km=header(3)
         close(IUNIT)
      else
         ncio = MAPL_NCIOOpen("fvcore_internal_restart_in")
         call MAPL_NCIOGetDimSizes(NCIO,lon=im,lat=jm,lev=km)
         call MAPL_NCIOClose(NCIO,destroy=.true.)
      end if 
   else
      call mpp_error(FATAL, 'ABORT: fvcore_internal_restart_in does not exist')
   endif

   if( file_exist("moist_internal_restart_in") ) then
      call MAPL_NCIOGetFileType("moist_internal_restart_in",ftype)
      if (ftype == 0) then
         isBinMoist = .false.
         NCIO = MAPL_NCIOOpen("moist_internal_restart_in")
         call MAPL_NCIOGetDimSizes(ncio,slices=nmoist,nVars=nVars)
         lcnt_var=2
         do ivar=1,nVars
            call MAPL_NCIOGetVarName(NCIO,ivar,vname)
            call MAPL_NCIOVarGetDims(NCIO,vname,ndims,dimSize)
            if (ndims==2) nmoist=nmoist-1
            if (ndims==3) then
            call moist_variables%push_back(trim(vname))
               if (trim(vname)=='Q') then
                  iq0=1
                  call tracer_names%insert(trim(vname),iq0)
               else
                  iq0=lcnt_var
                  call tracer_names%insert(trim(vname),iq0)
                  lcnt_var=lcnt_var+1
               end if
            end if
         enddo
         call MAPL_NCIOClose(ncio,destroy=.true.)
      else
         call rs_count( "moist_internal_restart_in",nmoist )
         if(mod(nmoist,km)/=0) then
           call mpp_error(FATAL, 'ABORT: '//'binary moist restart must have only 3D variabels')
         end if
         nVars = nmoist/km
         do ivar=1,nVars
            if (ivar < 10) then
               call tracer_names%insert(trim(moist_order(ivar)),ivar)
            else
               call tracer_names%insert("unknown",ivar)
            end if
         enddo
      end if
   else
      call mpp_error(FATAL, 'ABORT: moist_internal_restart_in does not exist')
   endif
   call mpp_broadcast(nmoist, mpp_root_pe())
   call mpp_broadcast(isBinMoist, mpp_root_pe())

   allocate ( r8_ak(npz+1) )
   allocate ( r8_bk(npz+1) )
   call set_eta(npz,ks,ptop,pint,r8_ak,r8_bk)
   FV_Atm(1)%ak = r8_ak
   FV_Atm(1)%bk = r8_bk
   deallocate ( r8_ak,r8_bk )
   nq = nmoist
   FV_Atm(1)%ncnst = nq/km
   if( is_master() ) then
      print *
      write(6,100)
100      format(2x,' k ','      A(k)    ',2x,' B(k)   ',2x,'  Pref    ',2x,'  DelP',/, &
            1x,'----',3x,'----------',2x,'--------',2x,'----------',2x,'---------' )
      k=1
      write(6,101) k,FV_Atm(1)%ak(k)*0.01, FV_Atm(1)%bk(k), FV_Atm(1)%ak(k)*0.01 + 1000.0*FV_Atm(1)%bk(k)
      do k=2,ubound(FV_Atm(1)%ak,1)
         write(6,102) k,FV_Atm(1)%ak(k)*0.01, FV_Atm(1)%bk(k), FV_Atm(1)%ak(k)*0.01 + 1000.0*FV_Atm(1)%bk(k), &
               (FV_Atm(1)%ak(k)-FV_Atm(1)%ak(k-1))*0.01 + 1000.0*(FV_Atm(1)%bk(k)-FV_Atm(1)%bk(k-1))
      enddo
      print *
101      format(2x,i3,2x,f10.6,2x,f8.4,2x,f10.4)
102      format(2x,i3,2x,f10.6,2x,f8.4,2x,f10.4,3x,f8.4)
103      format(2x,a,i6,3x,a,f7.2,a)
      write(6,103) 'Total Number of Tracers in  MOIST: ',nmoist ,'(/KM = ',float(nmoist) /float(km),')'
      print *
   endif

   do i=1,n_files

      if (file_exist(trim(extra_files(i)))) then
         call MAPL_NCIOGetFileType(trim(extra_files(i)),ftype)
         if (ftype ==0) then
            rst_files(i)%isBin=.false.
            rst_files(i)%file_name=trim(extra_files(i))
            NCIO = MAPL_NCIOOpen(trim(extra_files(i)))
            allocate(rst_files(i)%vars(size(ncio%vars)))
            edgeid=-1
            levid=-1
            call MAPL_NCIOGetDimSizes(ncio,lev=nLevs,edges=nEdges,levid=levid,edgeid=edgeid)
            do n=1,size(ncio%vars)
               rst_files(i)%vars(n)%name=trim(ncio%vars(n)%name)
               if (ncio%vars(n)%ndims==3) then
                  do j=1,size(ncio%vars(n)%dimids)
                     if (ncio%vars(n)%dimids(j) == edgeid) then
                        rst_files(i)%vars(n)%nlev= npz+1
                     else if (ncio%vars(n)%dimids(j) == levid) then
                        rst_files(i)%vars(n)%nlev= npz
                     end if
                  enddo
               else if (ncio%vars(n)%ndims==2) then
                  rst_files(i)%vars(n)%nlev=1
               end if
            enddo
            rst_files(i)%have_descriptor=.true.
            call MAPL_NCIOClose(NCIO,destroy=.true.)
         else
           call rs_count(trim(extra_files(i)),nlevs)
           if (mod(nlevs,km) /= 0) then
              rst_files(i)%have_descriptor=.false.
              nvars = 1
           else
              rst_files(i)%have_descriptor=.true.
              nvars=nlevs/km
           end if
           allocate(rst_files(i)%vars(nvars))
           rst_files(i)%file_name=trim(extra_files(i))
           rst_files(i)%isBin=.true.
           if (rst_files(i)%have_descriptor) then
              rst_files(i)%vars%nLev=npz
           else
              rst_files(i)%vars%nLev=nlevs
           end if
         end if
      else
         call mpp_error(FATAL, 'ABORT: '//trim(extra_files(i))//' does not exist')
      end if

   end do

   call print_memuse_stats('interp_restarts: begining get_external_ic')

   FV_Atm(1)%flagstruct%Make_NH = .false. ! Do this after rescaling
   if (jm == 6*im) then
      call get_geos_ic( FV_Atm, rst_files, .true.)
   else
      call get_geos_ic( FV_Atm, rst_files, .false.)
   endif
   FV_Atm(1)%flagstruct%Make_NH = .true. ! Reset this for later

   is = FV_Atm(1)%bd%isc
   ie = FV_Atm(1)%bd%iec
   js = FV_Atm(1)%bd%jsc
   je = FV_Atm(1)%bd%jec

   if (scale_rst) then
      call scale_drymass(fv_atm,tracer_names,rc=status)
      VERIFY_(status)
   end if

   if (FV_Atm(1)%flagstruct%Make_NH) then
      if (is_master()) print*, 'Updating FV3 NonHydrostatic State'
      do k=1,npz
      do j=js,je
      do i=is,ie
      FV_Atm(1)%w(i,j,k) = 0.0
      FV_Atm(1)%delz(i,j,k) = (-MAPL_RGAS/MAPL_GRAV)*FV_Atm(1)%pt(i,j,k)*(log(FV_Atm(1)%pe(i,k+1,j))-log(FV_Atm(1)%pe(i,k,j)))
      FV_Atm(1)%pkz(i,j,k)  = exp( MAPL_KAPPA*log((-MAPL_RGAS/MAPL_GRAV)*FV_Atm(1)%delp(i,j,k)*FV_Atm(1)%pt(i,j,k)*    &
                                       (1.0+(MAPL_RVAP/MAPL_RGAS - 1.)*FV_Atm(1)%q(i,j,k,1))/FV_Atm(1)%delz(i,j,k)) )
      enddo
      enddo
      enddo
   endif

   allocate(pt_local(is:ie,js:je,npz))
   pt_local=0.0d0
   do k=1,npz
! Convert to Potential Temperature
      pt_local(is:ie,js:je,k) = FV_Atm(1)%pt(is:ie,js:je,k)/FV_Atm(1)%pkz(is:ie,js:je,k)
   enddo



   npes_x=fv_atm(1)%layout(1)
   npes_y=fv_atm(1)%layout(2)
   isl=is
   iel=ie
   jsl=(npx-1)*(tile-1)+js
   jel=(npx-1)*(tile-1)+je

   call ArrDescrInit(Arrdes,MPI_COMM_WORLD,npx-1,(npx-1)*6,npes_x,npes_y*6,n_readers,n_writers,isl,iel,jsl,jel,rc=status)
   call ArrDescrSet(arrdes,offset=0_MPI_OFFSET_KIND)
   amWriter = arrdes%writers_comm/=MPI_COMM_NULL
   grid = simple_cs_grid_creator(npx-1,(npx-1)*6,npes_x,npes_y*6,status)
   call MPI_Info_create(info,status)

   call print_memuse_stats('interp_restarts: going to write restarts')

! write fvcore_internal_rst
   if( file_exist("fvcore_internal_restart_in") ) then

         write(fname1, "('fvcore_internal_rst_c',i4.4,'_',i3.3,'L')") npx-1,npz
         if (is_master()) print*, 'Writing : ', TRIM(fname1)
         if (isBinFV) then

            open(IUNIT,file='fvcore_internal_restart_in' ,access='sequential',form='unformatted',status='old')
            !if (AmWriter) open(OUNIT,file=TRIM(fname1),access='sequential',form='unformatted')
            if (n_writers==1) then
               OUNIT=getfile(TRIM(fname1),form='unformatted',rc=status)
               VERIFY_(status)
            else
               if (AmWriter) then
                  call MPI_FILE_OPEN(arrdes%writers_comm, fname1, MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                                      info, OUNIT, STATUS)
                  VERIFY_(STATUS)
               end if
            end if

            ! Headers
            read (IUNIT, IOSTAT=status) header 
            if(n_writers > 1) then
               call Write_Parallel(HEADER, OUNIT, ARRDES=ARRDES, RC=status)
               VERIFY_(STATUS)
            else
               if (amwriter) write(OUNIT) header
            endif
            if (is_master()) print*, header

            read (IUNIT, IOSTAT=status) header(1:5)
            if (is_master()) print*, header(1:5)  
            header(1) = (npx-1)
            header(2) = (npy-1)*6
            header(3) = npz

            if(n_writers > 1) then
               call Write_Parallel(HEADER(1:5), OUNIT, ARRDES=ARRDES, RC=status)
               VERIFY_(STATUS)
            else
               if (amwriter) write(OUNIT) header(1:5)
            endif

            if (is_master()) print*, header(1:5) 
            close(IUNIT)

         else

            ncio = MAPL_NCIOOpen("fvcore_internal_restart_in",rc=status)
            call MAPL_NCIOGetDimSizes(ncio,nVars=nVars)
            if (AmWriter) then
               imc = npx-1
               jmc = imc*6
               call MAPL_NCIOChangeRes(ncio,ncioOut,lonSize=imc,latSize=jmc,levSize=npz)
               if (.not.fv_atm(1)%flagstruct%hydrostatic) then
                  nfv_vars = 9
               else
                  nfv_vars = 7
               end if
               call MAPL_NCIOSet(ncioOut,filename=fname1,nVars=nfv_vars,overwriteVars=.true.)
               !start creating file. Can not simply use input because if it is lat-lon will not have dz or w
               call MAPL_NCIOGetDimSizes(ncioOut,lonid=lonid,latid=latid,levid=levid,edgeid=edgeid)
               call MAPL_NCIOAddVar(ncioOut,"AK",(/edgeid/),6,units="Pa",long_name="hybrid_sigma_pressure_a",rc=status)
               VERIFY_(status)
               call MAPL_NCIOAddVar(ncioOut,"BK",(/edgeid/),6,units="Pa",long_name="hybrid_sigma_pressure_b",rc=status)
               VERIFY_(status)
               call MAPL_NCIOAddVar(ncioOut,"U",(/lonid,latid,levid/),6,units="m s-1",long_name="eastward_wind",rc=status)
               VERIFY_(status)
               call MAPL_NCIOAddVar(ncioOut,"V",(/lonid,latid,levid/),6,units="m s-1",long_name="northward_wind",rc=status)
               VERIFY_(status)
               call MAPL_NCIOAddVar(ncioOut,"PT",(/lonid,latid,levid/),6,units="K Pa$^{-\kappa}$",long_name="scaled_potential_temperature",rc=status)
               VERIFY_(status)
               call MAPL_NCIOAddVar(ncioOut,"PE",(/lonid,latid,edgeid/),6,units="Pa",long_name="air_pressure",rc=status)
               VERIFY_(status)
               call MAPL_NCIOAddVar(ncioOut,"PKZ",(/lonid,latid,levid/),6,units="Pa$^\kappa$",long_name="pressure_to_kappa",rc=status)
               VERIFY_(status)
               
               ! if dz and w were not in the original file add them
               ! they need to be in there for the restart

               if (.not.fv_atm(1)%flagstruct%hydrostatic) then
                  call MAPL_NCIOAddVar(ncioOut,"DZ",(/lonid,latid,levid/),6,units="m",long_name="height_thickness",rc=status)
                  VERIFY_(status)
                  call MAPL_NCIOAddVar(ncioOut,"W",(/lonid,latid,levid/),6,units="m s-1",long_name="vertical_velocity",rc=status)
                  VERIFY_(status)
               endif
               call MAPL_NCIOCreateFile(ncioOut,comm=arrdes%writers_comm,info=info,rc=status)
               VERIFY_(status)
            end if

         end if

! AK and BK
         allocate ( r8_akbk(npz+1) )
         r8_akbk = FV_Atm(1)%ak
         if (isBinFV) then
            if (n_writers == 1) then
               if (AmWriter) write(OUNIT) r8_akbk
            else
               call write_parallel(r8_akbk,ounit,arrdes=arrdes,rc=status)
               VERIFY_(status)
            end if
         else  
            if (AmWriter) call MAPL_VarWrite(ncioOut,"AK",r8_akbk)
         end if
         r8_akbk = FV_Atm(1)%bk
         if (isBinFV) then
            if (n_writers == 1) then
               if (AmWriter) write(OUNIT) r8_akbk
            else
               call write_parallel(r8_akbk,ounit,arrdes=arrdes,rc=status)
               VERIFY_(status)
            end if
         else  
            if (AmWriter) call MAPL_VarWrite(ncioOut,"BK",r8_akbk)
         end if
         deallocate ( r8_akbk )

         allocate(r8_local(is:ie,js:je,npz+1))

! U
         if (is_master()) print*, 'Writing : ', TRIM(fname1), ' U'
         r8_local(is:ie,js:je,1:npz) = FV_Atm(1)%u(is:ie,js:je,1:npz)
         if (isBinFV) then
            if (n_writers==1) then
               call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz),rc=status)
               VERIFY_(status)
            else
               call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
               VERIFY_(status)
            end if
         else
            call MAPL_VarWrite(ncioOut,"U",r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
            VERIFY_(status)
         end if   
! V
         if (is_master()) print*, 'Writing : ', TRIM(fname1), ' V'
         r8_local(is:ie,js:je,1:npz) = FV_Atm(1)%v(is:ie,js:je,1:npz)
         if (isBinFV) then
            if (n_writers==1) then
               call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz),rc=status)
               VERIFY_(status)
            else
               call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
               VERIFY_(status)
            end if
         else
            call MAPL_VarWrite(ncioOut,"V",r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
               VERIFY_(status)
         end if   
! PT
         if (is_master()) print*, 'Writing : ', TRIM(fname1), ' PT'

         if (isBinFV) then
            if (n_writers==1) then
               call MAPL_VarWrite(OUNIT,grid,pt_local(is:ie,js:je,1:npz),rc=status)
               VERIFY_(status)
            else
               call MAPL_VarWrite(OUNIT,grid,pt_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
               VERIFY_(status)
            end if
         else
            call MAPL_VarWrite(ncioOut,"PT",pt_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
            VERIFY_(status)
         end if

! PE
         if (is_master()) print*, 'Writing : ', TRIM(fname1), ' PE'
         do k=1,npz+1
            r8_local(is:ie,js:je,k) = FV_Atm(1)%pe(is:ie,k,js:je)
         enddo
         if (isBinFV) then
            if (n_writers==1) then
               call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz+1),rc=status)
               VERIFY_(status)
            else
               call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz+1),arrdes=arrdes,rc=status)
               VERIFY_(status)
            end if
         else
            call MAPL_VarWrite(ncioOut,"PE",r8_local(is:ie,js:je,1:npz+1),arrdes=arrdes,rc=status)
            VERIFY_(status)
         end if
! PKZ
         if (is_master()) print*, 'Writing : ', TRIM(fname1), ' PKZ'
         r8_local(is:ie,js:je,1:npz) = FV_Atm(1)%pkz(is:ie,js:je,1:npz)
         if (isBinFV) then
            if (n_writers==1) then
               call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz),rc=status)
               VERIFY_(status)
            else
               call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
               VERIFY_(status)
            end if
         else
            call MAPL_VarWrite(ncioOut,"PKZ",r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
            VERIFY_(status)
         end if

         if (.not. fv_atm(1)%flagstruct%hydrostatic) then
! DZ
            if (is_master()) print*, 'Writing : ', TRIM(fname1), ' DZ'
            r8_local(is:ie,js:je,1:npz) = FV_Atm(1)%delz(is:ie,js:je,1:npz)
            if (isBinFV) then
               if (n_writers==1) then
                  call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz),rc=status)
                  VERIFY_(status)
               else
                  call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
                  VERIFY_(status)
               end if
            else
               call MAPL_VarWrite(ncioOut,"DZ",r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
               VERIFY_(status)
            end if

! W
            if (is_master()) print*, 'Writing : ', TRIM(fname1), ' W'
            if (is_master()) print*, 'Writing : ', TRIM(fname1), ' DZ'
            r8_local(is:ie,js:je,1:npz) = FV_Atm(1)%w(is:ie,js:je,1:npz)
            if (isBinFV) then
               if (n_writers==1) then
                  call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz),rc=status)
                  VERIFY_(status)
               else
                  call MAPL_VarWrite(OUNIT,grid,r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
                  VERIFY_(status)
               end if
            else
               call MAPL_VarWrite(ncioOut,"W",r8_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
               VERIFY_(status)
            end if
         endif

         if (isBinFV) then
            if (n_writers > 1) then
               if (AmWriter) then
                  call MPI_FILE_CLOSE(OUNIT,status)
                  VERIFY_(status)
               end if
            else
               close (OUNIT)
            end if
         else
            if (AmWriter) call MAPL_NCIOClose(ncioOut,destroy=.true.,rc=status)
         end if

         deallocate (r8_local)
         deallocate (pt_local)

      endif

! MOIST
!
      allocate(r4_local(is:ie,js:je,npz+1))
      allocate(r4_local2D(is:ie,js:je))

      if( file_exist("moist_internal_restart_in") ) then
         write(fname1, "('moist_internal_rst_c',i4.4,'_',i3.3,'L')") npx-1,npz
         if (is_master()) print*, 'Writing : ', TRIM(fname1)
         if (isBinMoist) then
            call ArrDescrSet(arrdes,offset=0_MPI_OFFSET_KIND)
            if (n_writers>1) then
               if (AmWriter) then
                  call MPI_FILE_OPEN(arrdes%writers_comm, fname1, MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                                      info, OUNIT, STATUS)
                  VERIFY_(STATUS)
               end if
            else
               ounit = getfile(trim(fname1),form='unformatted',rc=status)
               VERIFY_(status)
               !open(OUNIT,file=fname1,access='sequential',form='unformatted')
            end if
         else
            !if (AmWriter) then
               ncio = MAPL_NCIOOpen("moist_internal_restart_in",rc=status)
               imc = npx-1
               jmc = imc*6
               call MAPL_NCIOChangeRes(ncio,ncioOut,lonSize=imc,latSize=jmc,levSize=npz)
               call MAPL_NCIOSet(ncioOut,filename=fname1)
            if (AmWriter) then
               call MAPL_NCIOCreateFile(ncioOut,comm=arrdes%writers_comm,info=info,rc=status)
               VERIFY_(status)
               call MAPL_NCIOClose(ncio,destroy=.true.,rc=status)
               VERIFY_(status)
            end if
         end if
         lcnt_var=2
         ! binary path
         if (isBinMoist) then
            do iq=1,FV_Atm(1)%ncnst
               if (is_master()) print*, 'Writing : ', TRIM(fname1), ' ', iq
               r4_local(is:ie,js:je,1:npz) = FV_Atm(1)%q(is:ie,js:je,:,iq)
               iq0=iq
               r4_local(is:ie,js:je,1:npz) = FV_Atm(1)%q(is:ie,js:je,:,iq0)
               if (n_writers == 1) then
                  call MAPL_VarWrite(OUNIT,grid,r4_local(is:ie,js:je,1:npz),rc=status)
                  VERIFY_(status)
               else
                  call MAPL_VarWrite(OUNIT,grid,r4_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
                  VERIFY_(status)
               end if
            end do
            if (n_writers == 1) then
               close (OUNIT)
            else
               if (AmWriter) then
                  call MPI_FILE_CLOSE(OUNIT,status)
                  VERIFY_(status)
               end if
            end if
         else
            ! netcdf path
            call MAPL_NCIOGetDimSizes(ncioOut,nVars=nVars)
            lcnt_var=2
            do ivar=1,nVars
               call MAPL_NCIOGetVarName(ncioOut,ivar,vname)
               if (is_master()) print*, 'Writing : ', TRIM(fname1), ' ', ivar
               call MAPL_NCIOVarGetDims(ncioOut,vname,ndims,dimSize)
               if (ndims==2) then
                  r4_local2d(is:ie,js:je)=0.0
                  call MAPL_VarWrite(ncioOut,vname,r4_local2d(is:ie,js:je),arrdes=arrdes,rc=status)
                  VERIFY_(status)
               else if (ndims==3) then
                  if (trim(vname)=='Q') then
                     iq0=1
                  else
                     iq0=lcnt_var
                     lcnt_var=lcnt_var+1
                  end if
                  r4_local(is:ie,js:je,1:npz) = FV_Atm(1)%q(is:ie,js:je,:,iq0)
                  call MAPL_VarWrite(ncioOut,vname,r4_local(is:ie,js:je,1:npz),arrdes=arrdes,rc=status)
                  VERIFY_(status)
               end if
            end do
            if (AmWriter) call MAPL_NCIOClose(ncioOut,destroy=.true.,rc=status)
            !if (AmWriter) call MAPL_NCIOClose(ncioOut,destroy=.true.,rc=status)
         end if
      end if
      deallocate(r4_local)
 
! extra restarts
!
      do ifile=1,size(rst_files)

         if (is_master()) write(*,*)'Writing results of ',trim(rst_files(ifile)%file_name) 
         !nq = index(rst_files(ifile)%file_name,"restart")
         !fname2 = trim(rst_files(ifile)%file_name(1:nq-1))//"rst_c"
         !write(fname1, "(A,i4.4,'_',i3.3,'L')") trim(fname2),npx-1,npz
         fname1=extra_output(ifile)
         if (is_master()) print*, 'Writing : ', TRIM(fname1)
         call ArrDescrSet(arrdes,offset=0_MPI_OFFSET_KIND)
         if (rst_files(ifile)%isBin) then
            if (n_writers>1) then
               if (AmWriter) then
                  call MPI_FILE_OPEN(arrdes%writers_comm, fname1, MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                                      info, ounit, STATUS)
                  VERIFY_(STATUS)
               end if
            else
               ounit=getfile(trim(fname1),form='unformatted',rc=status)
               VERIFY_(status)
            end if
         else
            !if (AmWriter) then
               ncio = MAPL_NCIOOpen(trim(rst_files(ifile)%file_name),rc=status)
               imc = npx-1
               jmc = imc*6
               call MAPL_NCIOChangeRes(ncio,ncioOut,lonSize=imc,latSize=jmc,levSize=npz)
               call MAPL_NCIOSet(ncioOut,filename=fname1)
            if (AmWriter) then
               call MAPL_NCIOCreateFile(ncioOut,comm=arrdes%writers_comm,info=info,rc=status)
            end if
         end if
         do iq=1,size(rst_files(ifile)%vars)
            if (.not.rst_files(ifile)%isBin .and. AmWriter) then
               call MAPL_NCIOGetVarName(ncio,iq,vname) 
            end if
            nlev = rst_files(ifile)%vars(iq)%nlev
            allocate(r4_local(is:ie,js:je,nlev))
            if (rst_files(ifile)%isBin) then
               if (nlev/=1) then
                  r4_local(is:ie,js:je,1:nlev)=rst_files(ifile)%vars(iq)%ptr3d(is:ie,js:je,1:nlev)
                  if (n_writers == 1) then
                     call MAPL_VarWrite(ounit,grid,r4_local(is:ie,js:je,1:nlev),rc=status)
                     VERIFY_(status)
                  else
                     call MAPL_VarWrite(ounit,grid,r4_local(is:ie,js:je,1:nlev),arrdes=arrdes,rc=status)
                     VERIFY_(status)
                  end if
               else
                  r4_local2d(is:ie,js:je)=rst_files(ifile)%vars(iq)%ptr2d(is:ie,js:je)
                  if (n_writers == 1) then
                     call MAPL_VarWrite(ounit,grid,r4_local2d(is:ie,js:je),rc=status)
                     VERIFY_(status)
                  else
                     call MAPL_VarWrite(ounit,grid,r4_local2d(is:ie,js:je),arrdes=arrdes,rc=status)
                     VERIFY_(status)
                  end if
               end if
            else
               if (nlev/=1) then
         if (is_master()) print*, 'Writing 3D: ', iq, TRIM(vname)
                  r4_local(is:ie,js:je,1:nlev)=rst_files(ifile)%vars(iq)%ptr3d(is:ie,js:je,1:nlev)
                  call MAPL_VarWrite(ncioOut,vname,r4_local(is:ie,js:je,1:nlev),arrdes)
               else
         if (is_master()) print*, 'Writing 2D: ', TRIM(vname)
                  r4_local2d(is:ie,js:je)=rst_files(ifile)%vars(iq)%ptr2d(is:ie,js:je)
                  call MAPL_VarWrite(ncioOut,vname,r4_local2d(is:ie,js:je),arrdes=arrdes)
               end if
            end if
         end do
         deallocate(r4_local)
         if ( rst_files(ifile)%isBin) then 
           if (n_writers == 1) then
              close (ounit)
           else
              if (AmWriter) then
                 call MPI_FILE_CLOSE(ounit,status)
                 VERIFY_(status)
              end if
           end if
         else
           !if (AmWriter) then
              call MAPL_NCIOClose(ncio,destroy=.true.,rc=status)
           if (AmWriter) then
              call MAPL_NCIOClose(ncioOut,destroy=.true.,rc=status)
           end if
         end if 

      end do

      deallocate(r4_local2D)

! Finalize SHMEM in MAPL
   call MAPL_FinalizeShmem (rc=status)

   call fv_end(fv_atm, grids_on_this_pe, .false.)
   call fms_end()

contains

   subroutine rs_count( filename,nrecs )
      implicit none
      integer        nrecs,rc
      character(*)   filename

      open  (55,file=trim(filename),form='unformatted',access='sequential')
      nrecs =  0
      rc =  0
      do while (rc.eq.0)
         read (55,iostat=rc)
         if( rc.eq.0 ) nrecs = nrecs + 1
      enddo
      close (55)

      return
   end subroutine rs_count

end program interp_restarts
