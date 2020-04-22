PROGRAM loss_during_day
! ifort -c -openmp -r8 loss_during_day.f90
! ifort -c -openmp -r8 m_loss_during_routines.f90
! ifort -o loss_during_day -openmp -r8 m_loss_during_routines.o loss_during_day.o 
! !USAGE
!     loss_during_day -js jseg

  use loss_during_routines
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SM: 11-22-2013
! for very sandy classes (a_sand > 90) dtstep is set to 1.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT NONE
      character*128 :: usage = "loss_during_day -js jseg"
      integer, parameter :: nlay=500,IM=20,MAXITS=20,IOUT=50,NTYPS=253
      integer, parameter :: ndays=100
      real, parameter :: rzdep = 1.0,zthick=.01
      REAL, dimension (nlay)  ::  w,psi,ak0,ak,flux,flux0,weq,weqnw, &
           winit,delw,c0,c1,c2,c3,h,dhdw,dkdw,c1x,c2x,c3x,c0x
      integer, parameter :: islopeind = 51
      integer :: itab,ianom,iwtab,irz,i,k,itstep,nstot,ns,jseg,job
      real :: wtdep,wanom,z,wtot,wrzstart,wtotnw,wdif,bterm,bterm2
      real :: dmax,dopt,dopt2,error,dnew,term1,term2,wtest,errtest
      real :: depth,wrzeq,term3,term4,term5,rzanew,wlrz,frac,eps,wrexc
      real :: dtstep
      REAL :: SS(nlay,21)
      character*300 ofile,opath,path
      character*30 arg
      real, dimension(ntyps) ::  a_bee,a_psis,a_aksat,a_poros, &
           a_sand,a_clay,a_silt,a_oc
      character*9, dimension (ntyps) :: a_class
      real ::  bee,psis,aksat,poros,class,dummy,oc, atau,btau
      character*300 :: soilfile
      logical :: file_exists,mult_jobs
      real, dimension (81,12) :: rzw, sfexc,tscale

! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
!
! NOTE: "!$" is for conditional compilation
!
logical :: running_omp = .false.
!
!$ integer :: omp_get_thread_num, omp_get_num_threads
!
integer :: n_threads=1, li, ui, kt, n
integer, allocatable, dimension (:) :: unit_number
!
integer, dimension(:), allocatable :: low_ind, upp_ind
!
! ------------------------------------------------------------------
!data a_bee /3.30, 3.80, 4.34, 5.25, 3.63, 5.96, 7.32, 8.41, &
!     8.34, 9.70, 10.78, 12.93/
!data a_psis/-0.05, -0.07, -0.16, -0.65, -0.84, -0.24, -0.12, &
!     -0.63, -0.28, -0.12, -0.58, -0.27/
!data a_poros/0.373, 0.386, 0.419, 0.476, 0.471, 0.437, 0.412, &
!     0.478, 0.447, 0.415, 0.478, 0.450/
!data a_aksat/2.45e-05, 1.75e-05, 8.35e-06, 2.36e-06, 1.1e-06, &
!     4.66e-06, 6.31e-06, 1.44e-06, 2.72e-06, 4.25e-06, 1.02e-06, 1.33e-06/
!       
! ----------- OpenMP PARALLEL ENVIRONMENT ----------------------------
!
! FIND OUT WHETHER -omp FLAG HAS BEEN SET DURING COMPILATION
!
!$ running_omp = .true.         ! conditional compilation
!
! ECHO BASIC OMP VARIABLES
!
!$OMP PARALLEL DEFAULT(NONE) SHARED(running_omp,n_threads) 
!
!$OMP SINGLE
!
!$ n_threads = omp_get_num_threads()
!
!$ write (*,*) 'running_omp = ', running_omp
!$ write (*,*)
!$ write (*,*) 'parallel OpenMP with ', n_threads, 'threads'
!$ write (*,*)
!$OMP ENDSINGLE
!
!$OMP CRITICAL
!$ write (*,*) 'thread ', omp_get_thread_num(), ' alive'
!$OMP ENDCRITICAL
!
!$OMP BARRIER
!
!$OMP ENDPARALLEL

! Process Arguments
!------------------
    mult_jobs = .false.
    jseg = ntyps
    job=1
    I = iargc()

    if(I < 1 ) then
       print *, "Job Segment is not specified: ", i
       print *, trim(Usage)
       print *, 'Assumed a single job and continuing ....'
    end if

    if (i>0) then 
       call getarg(1,arg)
       if ( trim(arg) == '-js' ) then
          call getarg(2,arg)
          read(arg,'(i2)') job
          mult_jobs = .true.
          jseg =32
          print *,job,jseg
       endif
    endif

     allocate(unit_number (1:n_threads))
     unit_number = 10

     allocate(low_ind(n_threads))
     allocate(upp_ind(n_threads))
     low_ind(1)         = 1
     upp_ind(n_threads) = NTYPS
      
     if (running_omp)  then
        if(mult_jobs ) then
           low_ind(1) = (job-1)*jseg + 1
           upp_ind(n_threads) = job*jseg
           if(job == 8) upp_ind(n_threads) = NTYPS
        endif
        
        do i=1,n_threads-1
           upp_ind(i)   = low_ind(i) + (jseg/n_threads) - 1 
           low_ind(i+1) = upp_ind(i) + 1
           unit_number(i) = unit_number(i) + i  
        end do

        if(job == 8) then 

           low_ind (14) = 251
           low_ind (15) = 252
           low_ind (16) = 253

           upp_ind (14) = 251
           upp_ind (15) = 252
           upp_ind (16) = 253

        endif
     end if

     soilfile='/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/' &
           //'/Woesten_SoilParam/Soil_param_100_mineral_3_OC_026_046_112_Woesten_topsoil.txt'

     open (10, file=trim(soilfile),form='formatted',status='old', &
          action='read')
     read (10,'(a10)')

     do k=1,NTYPS
        read (10,*)i,a_class(k),a_sand(k),a_clay(k),a_silt(k),a_oc(k),dummy,    &
             a_poros(k),a_bee(k),a_psis(k),a_aksat(k)
!        print *,a_class(k),a_sand(k),a_clay(k),a_silt(k),a_oc(k),dummy,    &
!             a_poros(k),a_bee(k),a_psis(k),a_aksat(k)
     end do

     close (10,status='keep')

     path = '/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/' &
          //'/Woesten_SoilParam/loss_pd_top/'

!$OMP PARALLELDO DEFAULT(NONE)                                    &
!$OMP SHARED(A_BEE, A_PSIS,A_AKSAT,A_POROS,n_threads, low_ind,    &
!$OMP  upp_ind,unit_number,a_class,a_sand,a_clay,a_silt,a_oc,path)&
!$OMP PRIVATE(itab,ianom,iwtab,irz,i,k,itstep,nstot,ns,          &
!$OMP  wtdep,wanom,z,wtot,wrzstart,wtotnw,wdif,bterm,bterm2,     &
!$OMP  dmax,dopt,dopt2,error,dnew,term1,term2,wtest,errtest,     &
!$OMP  depth,wrzeq,term3,term4,term5,rzanew,wlrz,frac,eps,wrexc, &
!$OMP  SS,bee,psis,aksat,poros,kt,li,ui,n,w,psi,ak0,ak,flux,     &
!$OMP  flux0,weq,weqnw,winit,delw,c0,c1,c2,c3,h,dhdw,dkdw,c1x,   &
!$OMP  c2x,c3x,c0x,ofile,oc,file_exists,dtstep)

     do kt=1,n_threads

        li = low_ind(kt)
        ui = upp_ind(kt)

      do n=li,ui

      bee   = a_bee(n)
      psis  = a_psis(n)
      poros = a_poros(n)
      aksat = a_aksat(n)
      oc    = a_oc(n)
      dtstep = 5.
      if(a_sand(n) > 90.) dtstep = 1.
      write (ofile,'(a,i2.2,i2.2,i4.4)')'loss_perday_rz1m_',nint(a_sand(n)),nint(a_clay(n)),nint(100*a_oc(n))

      print *,unit_number(kt),trim(ofile)

      inquire(file=trim(path)//trim(ofile), exist=file_exists)

      if (file_exists) then
      
         print *,'File exists :',trim(ofile)
       
      else

      open(unit_number(kt),file=trim(path)//trim(ofile),form='formatted',status='unknown') 

      do itab=100,500,5
!         if(mod(itab,10).eq.0) write(*,*) 'Starting depth #',itab,kt
         wtdep=itab/100.
         do ianom=-20,20
            wanom=ianom/100.
            
! Start with equilibrium profile.  i=1 refers to lowest layer.

          iwtab=nlay-int(wtdep/zthick+.01)
          irz=nlay-int(rzdep/zthick+.01)

          do i=1,iwtab
             w(i)=1.
             winit(i)=1.
          enddo

          do i=iwtab+1,nlay
             z=(i-iwtab-.5)*zthick
             w(i)=((psis-z)/psis)**(-1./bee)
             winit(i)=((psis-z)/psis)**(-1./bee)
          enddo

          wtot=0.
          do i=1,nlay
             delw(i)=0.
             weq(i)=w(i)
             wtot=wtot+w(i)*poros*zthick
         enddo

! Impose anomaly in root zone
         wrzstart=0.
         do i=irz+1,nlay
            w(i)=w(i)+wanom
            if(w(i) .gt. 1.) w(i)=1.
            if(w(i) .lt. 0.) w(i)=0.
            wrzstart=wrzstart+w(i)*poros*zthick
         enddo
         
! Compute equilibrium profile, given anomaly:
         wtotnw=0.
         do i=1,nlay
            wtotnw=wtotnw+w(i)*poros*zthick
         enddo
         wdif=wtotnw-wtot
         
         bterm=bee/(bee-1.)
         bterm2=1./bterm
         
         dmax=(wtotnw-wtot)/poros
         dopt=-999.
         error=10000000.
         do k=-100,100
             dnew=wtdep-k*wtdep/100.
             term1=((psis-dnew)/psis)**bterm2
             term2=-bterm*poros*psis*(term1-1.)
             wtest=poros*(nlay*zthick-dnew)+term2
             errtest=abs(wtest-wtotnw)
             if(errtest.lt.error) then
                error=errtest
                dopt=dnew
             endif
          enddo
          
          dopt2=dopt
          do k=-100,100
             dnew=dopt+k*wtdep/10000.
             term1=((psis-dnew)/psis)**bterm2
             term2=-bterm*poros*psis*(term1-1.)
             wtest=poros*(nlay*zthick-dnew)+term2
             errtest=abs(wtest-wtotnw)
             if(errtest.lt.error) then
                error=errtest
                dopt2=dnew
             endif
          enddo

          dopt=dopt2
          do k=-100,100
             dnew=dopt+k*wtdep/1000000.
             term1=((psis-dnew)/psis)**bterm2
             term2=-bterm*poros*psis*(term1-1.)
             wtest=poros*(nlay*zthick-dnew)+term2
             errtest=abs(wtest-wtotnw)
             if(errtest.lt.error) then
                error=errtest
                dopt2=dnew
             endif
          enddo
          
          depth=dopt2
          
          if(depth.gt.rzdep) then
             term1=((psis-depth)/psis)**bterm2
             term2=((psis-depth+rzdep)/psis)**bterm2
             wrzeq=poros*psis*(term1-term2)*bterm*(-1.)
          else
             term1=((psis-depth)/psis)**bterm2
             wrzeq=-bterm*poros*psis*(term1-1)
             wrzeq=wrzeq+poros*(rzdep-depth)
          endif
          
          
! Compute water mass that must leave root zone to attain equilibrium
          wrexc=wrzstart-wrzeq

! ---------------------------------------------------------

          itstep=int(dtstep+.001)
          if(mod(3600,itstep) .ne. 0) then
             write(*,*) 'Time step length of ',dtstep, &
                  ' does not evenly divide hour'
             stop
          endif
          nstot=int(24.*3.*1200/itstep)

          do ns=1,nstot

 !        Establish psi, K, w-equil values:
 
             do i=1,nlay
                if(w(i) .lt. 0.) print *,'w = 0 at i =',i,' nstep=',n,ns,bee,psis,aksat
                if(w(i) .lt. 1.e-20)w(i) = 1.e-20
                psi(i)=psis*(w(i)**(-bee))
                ak0(i)=aksat*(w(i)**(2.*bee+3))
                h(i)=(i-.5)*zthick+psi(i)
             enddo
             do i=2,nlay
!            ak(i)=ak0(i)*ak0(i-1)/(ak0(i)+ak0(i-1))
                ak(i)=(ak0(i)+ak0(i-1))/2.
             enddo

!       Compute fluxes.  Flux is defined as positive downward:
            do i=2,nlay
               flux(i)=dtstep*ak(i)*(h(i)-h(i-1))/zthick
            enddo

!       Load arrays of matrix elements:

            do i=1,nlay
               dkdw(i)=0.5*aksat*(2*bee+3)*w(i)**(2*bee+2)
               dhdw(i)=psis*(-bee)*w(i)**(-bee-1.)
            enddo

            do i=2,nlay-1
               c0(i)=(flux(i)-flux(i+1))/(poros*zthick)
               term1=dtstep/(zthick*zthick*poros)
               term2=(h(i+1)-h(i))*dkdw(i+1)
               term3=((ak0(i+1)+ak0(i))/2.)*dhdw(i+1)
               c1(i)=term1*(term2+term3)
               term2=(h(i+1)-h(i))*dkdw(i)
               term3=((ak0(i)+ak0(i+1))/2.)*dhdw(i)
               term4=(h(i)-h(i-1))*dkdw(i)
               term5=((ak0(i)+ak0(i-1))/2.)*dhdw(i)
               c2(i)=(term1*(term2-term3))-(term1*(term4+term5))-1.
               term2=(h(i)-h(i-1))*dkdw(i-1)
               term3=((ak0(i-1)+ak0(i))/2.)*dhdw(i-1)
               c3(i)=-term1*(term2-term3)
            enddo

            c0(nlay)=flux(nlay)/(poros*zthick)
            c1(nlay)=0.
            term1=dtstep/(zthick*zthick*poros)
            term2=(h(nlay)-h(nlay-1))*dkdw(nlay)
            term3=((ak0(nlay-1)+ak0(nlay))/2.)*dhdw(nlay)
            c2(nlay)=-1.-term1*(term2+term3)
            term2=(h(nlay)-h(nlay-1))*dkdw(nlay-1)
            term3=((ak0(nlay-1)+ak0(nlay))/2.)*dhdw(nlay-1)
            c3(nlay)=-term1*(term2-term3)
            
            c0(1)=flux(2)/(poros*zthick)
            term1=dtstep/(zthick*zthick*poros)
            term2=(h(2)-h(1))*dkdw(2)
            term3=((ak0(2)+ak0(1))/2.)*dhdw(2)
            c1(1)=-term1*(term2+term3)
            term2=(h(2)-h(1))*dkdw(1)
            term3=((ak0(2)+ak0(1))/2.)*dhdw(1)
            c2(1)=1.-(term1*(term2-term3))
            c3(1)=0.
            
            if(ns.eq.1) then
               do i=1,nlay
                  c0x(i)=c0(i)
                  c1x(i)=c1(i)
                  c2x(i)=c2(i)
                  c3x(i)=c3(i)
               enddo
               call tridag(c1x,c2x,c3x,c0x,delw,nlay)
            endif
            
            eps=.00001
            call GMRES (nlay, IM, c0, delw, SS, eps, MAXITS, IOUT, &
                 c1,c2,c3)
            

!       Update reservoirs
            do i=1,nlay
               w(i)=w(i)+delw(i)
               if(w(i) .gt.1.) then
                  delw(i+1)=delw(i+1)+(w(i)-1.)
                  w(i)=1.
               endif
            enddo
!           write(15,1015) w
!           format(1x,8e16.8)

            enddo  ! time-steps per 24.*3.*1200s loop


!c          stop
!c**** Compute water (wlrz) that actually left root zone
          rzanew=0.
          do i=irz+1,nlay
             rzanew=rzanew+w(i)*zthick*poros
          enddo
          wlrz=wrzstart-rzanew
          
          frac=wlrz/(wrexc+1.e-20)
          if(abs(wrexc).lt.1.e-5) frac=1.
          write(unit_number(kt),1000) wtdep,wanom,wrexc,frac
1000      format(1x,4e16.8)
          
       enddo      ! Anomaly loop
    enddo        ! Depth loop
    close (unit_number(kt), status='keep')
 endif
end do
end do

!$OMP ENDPARALLELDO

stop
end PROGRAM loss_during_day

