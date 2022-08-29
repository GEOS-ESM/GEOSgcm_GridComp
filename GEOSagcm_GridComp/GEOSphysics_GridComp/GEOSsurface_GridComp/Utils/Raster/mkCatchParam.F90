PROGRAM mkCatchParam

! !INTERFACE:
!
! !ARGUMENTS:
!
!  Usage = "mkCatchParam -x nx -y ny -g Gridname -b DL -v LBSV -e EASE"       
!     -x: Size of longitude dimension of input raster. DEFAULT: 8640
!     -y: Size of latitude dimension of input raster.  DEFAULT: 4320
!     -b: position of the dateline in the first box. DEFAULT: DC 
!     -g: Gridname  (name of the .til or .rst file without file extension)  
!     -v: LBCSV : Choose bcs version (ICA, NL3, NL4, NL5, or development)             
!     -e: EASE : This is optional if catchment.def file is available already or                    
!         the til file format is pre-Fortuna-2.                                                    
!     
!
! This program is good to generate  
! model, vegetation, soil, and MODIS albedo parameter files for the 
! catchment model implementation
!  
! Sarith Mahanama - March 23, 2012 
! Email: sarith.p.mahanama@nasa.gov

  use rmTinyCatchParaMod
  use process_hres_data
  use write_CatchParamsMod

  implicit none
  
  character*4          :: LBSV = 'DEF'
  character*128        :: GridName = ''
  character*128        :: ARG, MaskFile
  character*256        :: CMD
  character*1          :: opt
  character*4          :: EASE ='    '
  character*2          :: DL ='DC'    
  integer              :: I, J, iargc, nxt
  character(len=400), dimension (8) ::  Usage 
  logical              :: ease_grid, file_exists
  integer :: nc, nr, status
  

  USAGE(1) ="Usage: mkCatchParam -x nx -y ny -g Gridname -b DL -v LBCSV                                       "
  USAGE(2) ="     -x: Size of longitude dimension of input raster. DEFAULT: 8640                                     "
  USAGE(3) ="     -y: Size of latitude dimension of input raster.  DEFAULT: 4320                                     "
  USAGE(4) ="     -g: Gridname  (name of the .til or .rst file without file extension)                               "
  USAGE(5) ="     -b: Position of the dateline in the first grid box (DC or DE). DEFAULT: DC                         "
  USAGE(6) ="     -v  LBCSV : Choose bcs version (F25, GM4, ICA, NL3, NL4, NL5, or DEV)                              "

! Process Arguments                            
!------------------ 

  CALL get_command (cmd)
  inquire(file='clsm/mkCatchParam.log', exist=file_exists)
  if(file_exists) then
     open (log_file, file ='clsm/mkCatchParam.log', status='old', position='append', form='formatted',action='write')
  else
     open (log_file, file ='clsm/mkCatchParam.log', status='unknown', form='formatted',action='write')
     write (log_file,'(a)')trim(cmd)
     write (log_file,'(a)')' '
  endif

  I = iargc()
  if(I < 1 .or. I > 10) then
     write (log_file,'(a)') "Wrong Number of arguments: ", i
     do j = 1,size(usage)
        print "(sp,a100)", Usage(j)
     end do
     call exit(1)
  end if

  ease_grid = .false.
  nxt = 1
  call getarg(nxt,arg)
  do while(arg(1:1)=='-')
     opt=arg(2:2)
     if(len(trim(arg))==2) then
        if(scan(opt,'zh')==0) then
           nxt = nxt + 1
           call getarg(nxt,arg)
        endif
     else
        arg = arg(3:)
     end if
     select case (opt)
     case ('x')
        read(arg,'(i6)') nc
     case ('y')
        read(arg,'(i6)') nr
     case ('g')
        GridName = trim(arg)
     case ('v')
        LBSV = trim(arg)
        call init_bcs_config (trim(LBSV))
     case ('b')
        DL = trim(arg)
     case ('e')
        EASE = trim(arg)
        if(EASE=='EASE') ease_grid=.true.
     case default
        do j = 1,size(usage)
           print "(sp,a100)", Usage(j)
        end do
        call exit(1)
     end select
     nxt = nxt + 1
     call getarg(nxt,arg)
  end do

  call getenv ("MASKFILE"        ,MaskFile        )
 
  if(trim(Gridname) == '') then
      write (log_file,'(a)')'Unable to create parameters without til/rst files.... !'
      stop
   endif

   call write_CatchParams(nc, nr, GridName, DL, LBSV, MaskFile, rc = status)

END PROGRAM mkCatchParam
