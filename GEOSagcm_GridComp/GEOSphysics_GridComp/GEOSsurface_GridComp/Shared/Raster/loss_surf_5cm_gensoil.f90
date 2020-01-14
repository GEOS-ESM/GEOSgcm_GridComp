PROGRAM loss_surf_5cm_gensoil
! ifort -c -openmp -r8 m_loss_during_routines.f90
! ifort -c -openmp -r8 loss_surf_5cm_gensoil.f90
! ifort -o loss_surf_5cm_gensoil -openmp -r8 m_loss_during_routines.o loss_surf_5cm_gensoil.o 

  use loss_during_routines

      IMPLICIT NONE
      integer, parameter :: nlay =100,IM=20,MAXITS=20,IOUT=50,NTYPS=253
      real, parameter :: rzdep = 1.0,zthick=.01,dtstep=1.
      integer, parameter :: islopeind = 51
      real, dimension (nlay) ::  w,psi,ak0,ak,flux,flux0,weq,weqnw,  &
           winit,delw,c0,c1,c2,c3,h,dhdw,dkdw,c1x,c2x,c3x,c0x,wchk
      integer :: itab,ianom,iwtab,irz,i,k,itstep,nstot,ns
      real :: rzbot,psi0,zmz0,wtot,wrzstart,wtotnw,wdif,bterm,bterm2
      real :: dmax,dopt,error,dnew,term1,wtest,errtest,wrexc
      real :: wtotchk,wanom,term2,w1eq,w2eq,w3eq,w4eq,w5eq
      real :: term3,term4,term5,rzanew,wlrz,frac,eps
      real :: SS(nlay,21)
      character*300 sfile,opath,path
      real, dimension(NTYPS) ::  a_bee,a_psis,a_aksat,a_poros, &
           a_sand,a_clay,a_silt,a_oc,a_wp,a_wpsurf,a_porosurf
      character*9, dimension (ntyps) :: a_class
      real ::  bee,psis,aksat,poros,dummy,oc, atau,btau,atau_2cm,btau_2cm
      character*300 :: soilfile
      real, dimension (81,12) :: rzw, sfexc,tscale
      logical :: skip = .false.
      logical :: layer_2cm = .false.

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
!      data a_bee /3.30, 3.80, 4.34, 5.25, 3.63, 5.96, 7.32, 8.41, &
!           8.34, 9.70, 10.78, 12.93/
!      data a_psis/-0.05, -0.07, -0.16, -0.65, -0.84, -0.24, -0.12, &
!           -0.63, -0.28, -0.12, -0.58, -0.27/
!      data a_poros/0.373, 0.386, 0.419, 0.476, 0.471, 0.437, 0.412, &
!           0.478, 0.447, 0.415, 0.478, 0.450/
!      data a_aksat/2.45e-05, 1.75e-05, 8.35e-06, 2.36e-06, 1.1e-06, &
!           4.66e-06, 6.31e-06, 1.44e-06, 2.72e-06, 4.25e-06, 1.02e-06, 1.33e-06/

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

     allocate(unit_number (1:n_threads))
     unit_number = 10
     
!     -----------------------------------------------
!     Loop over root zone amounts and surface anomalies
!      if(n < 10) 
     allocate(low_ind(n_threads))
     allocate(upp_ind(n_threads))
     low_ind(1)         = 1
     upp_ind(n_threads) = NTYPS
      
     if (running_omp)  then
      do i=1,n_threads-1
        
         upp_ind(i)   = low_ind(i) + (NTYPS/n_threads) - 1 
         low_ind(i+1) = upp_ind(i) + 1
         unit_number(i) = unit_number(i) + i
      end do 
     end if

     soilfile='/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/' &
           //'/Woesten_SoilParam/Soil_param_100_mineral_3_OC_026_046_112_Woesten_topsoil.txt'

     open (10, file=trim(soilfile),form='formatted',status='old', &
          action='read')
     read (10,'(a10)')
 
     do k=1,NTYPS
        read (10,*)i,a_class(k),a_sand(k),a_clay(k),a_silt(k),a_oc(k),a_wpsurf(k),    &
             a_poros(k),a_bee(k),a_psis(k),a_aksat(k)
        a_porosurf(k) = a_poros(k)
      end do
     close (10,status='keep')

      
     path = '/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/' &
          //'//Woesten_SoilParam/loss_ph/'

if(.not.skip) then
!$OMP PARALLELDO DEFAULT(NONE)                                    &
!$OMP SHARED(A_BEE, A_PSIS,A_AKSAT,A_POROS,n_threads, low_ind,    &
!$OMP  upp_ind,unit_number,a_class,a_sand,a_clay,a_silt,a_oc,path &
!$OMP  ,layer_2cm)&
!$OMP PRIVATE(itab,ianom,iwtab,irz,i,k,itstep,nstot,ns,          &
!$OMP  wanom,wtot,wrzstart,wtotnw,wdif,bterm,bterm2,             &
!$OMP  dmax,dopt,error,dnew,term1,term2,wtest,errtest,           &
!$OMP  term3,term4,term5,rzanew,wlrz,frac,eps,wrexc,             &
!$OMP  SS,bee,psis,aksat,poros,kt,li,ui,n,w,psi,ak0,ak,flux,     &
!$OMP  flux0,weq,weqnw,winit,delw,c0,c1,c2,c3,h,dhdw,dkdw,c1x,   &
!$OMP  c2x,c3x,c0x,opath, rzbot,psi0,zmz0,w1eq,w2eq,w3eq,w4eq,   &
!$OMP  w5eq,wtotchk,wchk,oc)

     do kt=1,n_threads

        li = low_ind(kt)
        ui = upp_ind(kt)

      do n=li,ui

      bee   = a_bee(n)
      psis  = a_psis(n)
      poros = a_poros(n)
      aksat = a_aksat(n)
      oc    = a_oc(n)

      write (opath,'(a,i2.2,i2.2,i4.4)')'loss_perhour_surf_5cm_',nint(a_sand(n)),nint(a_clay(n)),nint(100*a_oc(n))
      if(layer_2cm) write (opath,'(a,i2.2,i2.2,i4.4)')'loss_perhour_surf_2cm_',nint(a_sand(n)),nint(a_clay(n)),nint(100*a_oc(n))
      open(unit_number(kt),file=trim(path)//trim(opath),form='formatted',status='unknown')


      do itab=20,100

        if(mod(itab,10).eq.0) write(*,*) 'Starting depth #',itab
        rzbot=itab/100.

        do ianom=-22,22,4
          wanom=ianom/40.

! Start with equilibrium profile.  i=1 refers to lowest layer.

          w(1)=rzbot
          psi0=psis*(w(1)**(-bee))
          do i=2,100
             zmz0=(i-1.)*zthick
             w(i)=((psi0-zmz0)/psis)**(-1./bee)
             winit(i)=w(i)
          enddo

          wtot=0.
          do i=1,nlay
            delw(i)=0.
            weq(i)=w(i)
            wtot=wtot+w(i)*poros*zthick
            enddo

!        Impose anomaly in surface layer
          wrzstart=0.
          k = 4
          if(layer_2cm) k = 1
          do i=nlay-k,nlay
            w(i)=w(i)+wanom
            if(w(i) .gt. 1.) w(i)=1.
            if(w(i) .lt. 0.001) w(i)=0.001
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

          dopt=-999.
          error=10000000.
          do k=-nlay,nlay
            dnew=psi0-k
            if(dnew.gt.-.001) dnew=-.001
            term1=((dnew-1.)/psis)**bterm2
            wtest=-bterm*poros*psis*1.*(term1-(dnew/psis)**bterm2)
            errtest=abs(wtest-wtotnw)
            if(errtest.lt.error) then
              error=errtest
              dopt=dnew
              endif
            enddo

          error=10000000.
          do k=-nlay,nlay
            dnew=dopt+k/100.
            if(dnew.gt.-.001) dnew=-.001
            term1=((dnew-1.)/psis)**bterm2
            wtest=-bterm*poros*psis*1.*(term1-(dnew/psis)**bterm2)
            errtest=abs(wtest-wtotnw)
            if(errtest.lt.error) then
              error=errtest
              dopt=dnew
              endif
            enddo

          error=10000000.
          do k=-nlay,nlay
            dnew=dopt+k/10000.
            if(dnew.gt.-.001) dnew=-.001
            term1=((dnew-1.)/psis)**bterm2
            wtest=-bterm*poros*psis*1.*(term1-(dnew/psis)**bterm2)
            errtest=abs(wtest-wtotnw)
            if(errtest.lt.error) then
              error=errtest
              dopt=dnew
              endif
            enddo

          error=10000000.
          do k=-nlay,nlay
            dnew=dopt+k/1000000.
            if(dnew.gt.-.001) dnew=-.001
            term1=((dnew-1.)/psis)**bterm2
            wtest=-bterm*poros*psis*1.*(term1-(dnew/psis)**bterm2)
            errtest=abs(wtest-wtotnw)
            if(errtest.lt.error) then
              error=errtest
              dopt=dnew
              endif
            enddo

          w1eq=((dopt-.995)/psis)**(-1./bee)
          w2eq=((dopt-.985)/psis)**(-1./bee)
          w3eq=((dopt-.975)/psis)**(-1./bee)
          w4eq=((dopt-.965)/psis)**(-1./bee)
          w5eq=((dopt-.955)/psis)**(-1./bee)

!**** Compute water mass that must leave root zone to attain equilibrium
          wrexc=wrzstart-(w1eq+w2eq+w3eq+w4eq+w5eq)*zthick*poros
          if(layer_2cm)wrexc=wrzstart-(w1eq+w2eq)*zthick*poros 
!**** check
          wtotchk=0.
          do i=1,100
            zmz0=(i-1.+.5)*zthick
            wchk(i)=((dopt-zmz0)/psis)**(-1./bee)
            wtotchk=wtotchk+wchk(i)*poros*zthick
         enddo

!**** ---------------------------------------------------------

          itstep=int(dtstep+.001)

          nstot=int(.001+3.*1200./dtstep)
!          nstot=1200/itstep

          do ns=1,nstot

!           Establish psi, K, w-equil values:
             do i=1,nlay
                if(w(i) .eq. 0.) print *,'w = 0 at i =',i,' nstep=',ns, bee,psis,poros,aksat
                psi(i)=psis*(w(i)**(-bee))
                ak0(i)=aksat*(w(i)**(2.*bee+3))
                h(i)=(i-.5)*zthick+psi(i)
             enddo

             do i=2,nlay
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
                if(w(i) .lt.0.001) then
                   delw(i+1)=delw(i+1)-(.001-w(i))
                   w(i)=0.001
                endif
             enddo
          enddo  ! time-steps per 24.*3.*1200s loop
           
          !          stop
          ! Compute water (wlrz) that actually left root zone
          rzanew=0.

          k = 4
          if(layer_2cm) k = 1

          do i=nlay-k,nlay
             rzanew=rzanew+w(i)*zthick*poros
          enddo
          wlrz=wrzstart-rzanew
          
          frac=wlrz/(wrexc+1.e-20)
          if(abs(wrexc).lt.1.e-5) frac=1.
          write(unit_number(kt),1000) wtot,wanom,wdif,frac
1000      format(1x,4e16.8)
          
       enddo      ! Anomaly loop
    enddo        ! Depth loop
    close (unit_number(kt), status='keep')
 end do
end do

!$OMP ENDPARALLELDO
if(layer_2cm) STOP
endif
   opath   = '/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/' & 
          //'//SoilClasses-SoilHyd-TauParam.dat'
   soilfile='/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/' &
        //'/Woesten_SoilParam/Soil_param_100_mineral_3_OC_026_046_112_Woesten_topsoil.txt'

   open (20, file=trim(soilfile),form='formatted',status='old', &
        action='read')
   read (20,'(a10)')

   open (11, file=trim(opath), form='formatted',status='unknown', &
        action = 'write')
   write (11,'(a)') 'Sand[w%] Clay[w%] Silt[w%] OC[w%] b[-] Psi_s[m] theta_s[m3/m3] wp[m3.m-3] K_s[ms-1] ATAU5 BTAU5 WP_SRF PORO_SRF ATAU_2cm BTAU_2cm' 



do n=1,NTYPS

! Computing atau_5cm and btau_5cm
   write (opath,'(a,i2.2,i2.2,i4.4)')'loss_perhour_surf_5cm_',nint(a_sand(n)),nint(a_clay(n)),nint(100*a_oc(n))
   open(10,file=trim(path)//trim(opath),form='formatted',status='old',action='read')
      do itab=1,81
        do ianom=1,12   
           read(10,1000) wtot,wanom,wdif,frac
           rzw  (itab,ianom) = wtot
           sfexc(itab,ianom) = wdif
           dummy             = -alog(1. - frac)
           tscale(itab,ianom)= 1./dummy
        end do
      end do
   close (10,status='keep')

   atau=(0.5*(tscale(islopeind,6)+tscale(islopeind,7)))*    &
        (0.5*(rzw(islopeind,6)+rzw(islopeind,7)))**3.

! estimate slope
   itstep=1
   do k=1,81 
      if(tscale(k,9) > 0.5*(tscale(islopeind,6)+tscale(islopeind,7))) itstep = k
   end do

   term3 = (0.5*(tscale(islopeind,6)+tscale(islopeind,7)) - tscale(itstep,9))  &
        /(tscale(itstep+1,9)-tscale(itstep,9))
   term4 = rzw(itstep,9)+term3*(rzw(itstep+1,9)-rzw(itstep,9))

   btau =-1.*(sfexc(islopeind,9)- 0.5*(sfexc(islopeind,6)+sfexc(islopeind,7)))   &
         /(term4 - 0.5*(rzw(islopeind,6)+rzw(islopeind,7)))

! Computing atau_2cm and btau_2cm
   write (opath,'(a,i2.2,i2.2,i4.4)')'loss_perhour_surf_2cm_',nint(a_sand(n)),nint(a_clay(n)),nint(100*a_oc(n))
   open(10,file=trim(path)//trim(opath),form='formatted',status='old',action='read')
      do itab=1,81
        do ianom=1,12   
           read(10,1000) wtot,wanom,wdif,frac
           rzw  (itab,ianom) = wtot
           sfexc(itab,ianom) = wdif
           dummy             = -alog(1. - frac)
           tscale(itab,ianom)= 1./dummy
        end do
      end do
   close (10,status='keep')

   atau_2cm=(0.5*(tscale(islopeind,6)+tscale(islopeind,7)))*    &
        (0.5*(rzw(islopeind,6)+rzw(islopeind,7)))**3.

! estimate slope
   itstep=1
   do k=1,81 
      if(tscale(k,9) > 0.5*(tscale(islopeind,6)+tscale(islopeind,7))) itstep = k
   end do

   term3 = (0.5*(tscale(islopeind,6)+tscale(islopeind,7)) - tscale(itstep,9))  &
        /(tscale(itstep+1,9)-tscale(itstep,9))
   term4 = rzw(itstep,9)+term3*(rzw(itstep+1,9)-rzw(itstep,9))

   btau_2cm =-1.*(sfexc(islopeind,9)- 0.5*(sfexc(islopeind,6)+sfexc(islopeind,7)))   &
         /(term4 - 0.5*(rzw(islopeind,6)+rzw(islopeind,7)))


   read (20,*)i,a_class(n),a_sand(n),a_clay(n),a_silt(n),a_oc(n),a_wp(n),    &
             a_poros(n),a_bee(n),a_psis(n),a_aksat(n)
   write (11,'(4f7.3,4f8.4,e13.5,2f12.7,2f8.4,4f12.7)')a_sand(n),a_clay(n),a_silt(n),a_oc(n),a_bee(n),a_psis(n), &
         a_poros(n),a_wp(n),a_aksat(n),atau,btau,a_wpsurf(n),a_porosurf(n),atau_2cm,btau_2cm

end do
close (11,status='keep')   
close (20,status='keep')
stop
end PROGRAM loss_surf_5cm_gensoil

