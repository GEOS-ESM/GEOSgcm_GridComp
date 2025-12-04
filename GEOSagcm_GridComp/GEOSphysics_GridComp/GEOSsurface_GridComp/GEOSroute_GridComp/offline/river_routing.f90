program main

use omp_lib        ! OpenMP library for parallel computing
use reservoir      ! Module for reservoir operations
use lake           ! Module for lake operations
use river_io       ! Module for river input/output

implicit none

! Define parameters and constants
real*8, parameter :: small = 1.D-48           ! A small value threshold for numerical comparisons
integer, parameter :: step_start = 9221       ! Start timestep (represents 1990-01-01)
integer, parameter :: step_end = 9226         ! End timestep (adjusted for different ranges)
logical :: is_coldstart = .True.              ! Logical flag for cold start
integer, parameter :: ny = 33                 ! Number of years (33 years)

real*8, parameter :: fac_kstr = 0.01D0       ! Factor for local stream scaling
real*8, parameter :: M = 0.45D0               ! Parameter in hydraulic geometry formula
real*8, parameter :: mm = 0.35D0              ! Parameter in hydraulic geometry formula

real*8, parameter :: dt = 86400.D0            ! Time step in seconds (1 day)

integer, parameter :: nmax = 373              ! Maximum number of catchments in a river
integer, parameter :: upmax = 34              ! Maximum number of upstream basins
integer, parameter :: nc = 291284             ! Total number of river cells
real*8, parameter :: rho = 1.D3               ! Water density in kg/m^3

! Declare variables
integer :: i, j, n, iter                      ! Loop indices and iteration variable

! Allocate dynamic arrays for variables
integer, allocatable, dimension(:) :: nts     ! Array for timestep indices
real*8, allocatable, dimension(:) :: qstr_clmt, qri_clmt, qin_clmt, & 
                                    llc_ori, llc, lstr, & 
                                    Qrunf, nume, deno, & 
                                    alp_s, alp_r, K, Kstr
real*8, allocatable, dimension(:) :: Ws, Wr   ! Water storage arrays for stream and river
real*8, allocatable, dimension(:) :: Qs0, ks, Ws_last, Qs, & 
                                    Qr0, kr, Cl, Al
real*8, allocatable, dimension(:) :: C1, C2, Qout, Qin, A1, P, B1, B2
integer, allocatable, dimension(:) :: nup     ! Number of upstream nodes
integer, allocatable, dimension(:,:) :: upID  ! IDs of upstream cells
real*8 :: co1, co2, co3                       ! Coefficients used in calculations
integer :: ui                                 ! Temporary upstream index variable

real*8, allocatable, dimension(:) :: lon, lat ! Longitude and latitude arrays

! Reservoir module variables
logical, parameter :: use_res = .True.        ! Flag to enable reservoir module
integer, parameter :: nres = 7250             ! Number of reservoirs
integer, allocatable, dimension(:) :: active_res, fld_res, cat2res  ! Reservoir attributes
real*8, allocatable, dimension(:) :: Wr_res, Q_res, cap_res, Qavg_res, ai_res, Qfld_thres, wid_res
integer, allocatable, dimension(:) :: type_res  ! Type of reservoir (0=inactive, 1-7=different functions)
real*8, allocatable, dimension(:,:) :: irr_sea_frac ! Irrigation and sea fraction for reservoirs

! Lake module variables
logical, parameter :: use_lake = .True.       ! Flag to enable lake module
integer, parameter :: nlake = 3917            ! Number of lakes
integer, allocatable, dimension(:) :: active_lake   ! Active lake flag
real*8, allocatable, dimension(:) :: area_lake, Wr_lake, Q_lake  ! Lake attributes

! Time-related variables
integer,dimension(ny) :: days_in_year=(/365,365,366,365,&
                                        365,365,366,365,&
                                        365,365,366,365,&
                                        365,365,366,365,&
                                        365,365,366,365,&
                                        365,365,366,365,&
                                        365,365,366,365,&
                                        365,365,366,365,365/) ! Number of days per year from 1990 to 2020
integer :: days_acc_year(ny), days_acc_noleap(12), days_acc_leap(12)  ! Accumulated days for leap and non-leap years
integer :: yr_cur, mon_cur, day_cur, d_res, step_prev  ! Current date variables and previous step
character(len=50) :: yr_s, mon_s, day_s               ! Year, month, day strings
character(len=500) :: inputdir                        ! Input directory path

! Allocate memory for variables
allocate(nts(nc))
allocate(qstr_clmt(nc), qri_clmt(nc), qin_clmt(nc))
allocate(llc_ori(nc), llc(nc), lstr(nc))
allocate(Qrunf(nc), nume(nc), deno(nc), alp_s(nc), alp_r(nc))
allocate(Ws(nc), Wr(nc))
allocate(Qs0(nc), ks(nc), Ws_last(nc), Qs(nc))
allocate(Qr0(nc), kr(nc), Cl(nc), Al(nc))
allocate(C1(nc), C2(nc), Qout(nc), Qin(nc), A1(nc), P(nc), B1(nc), B2(nc))
allocate(nup(nc))
allocate(upID(upmax,nc))
allocate(K(nc), Kstr(nc))

! Read input data
call read_input(nc, ny, upmax, days_in_year, fac_kstr, qstr_clmt, qri_clmt, nts, upID, nup, llc_ori, lstr, qin_clmt, K, Kstr, days_acc_year, days_acc_noleap, days_acc_leap, inputdir)

! Initialize reservoir module
call res_init(inputdir, nres, nc, use_res, active_res, Wr_res, Q_res, type_res, cap_res, Qavg_res, ai_res, fld_res, Qfld_thres, irr_sea_frac, cat2res, wid_res)

! Initialize lake module
call lake_init(inputdir, use_lake, nc, nlake, nres, active_res, active_lake, area_lake, Wr_lake, Q_lake)

! Calculate llc (length of river channel)
nume = qri_clmt**(2.D0-M) - qin_clmt**(2.D0-M)  ! Numerator for the llc calculation
deno = (2.D0-M) * (qri_clmt - qin_clmt) * (qri_clmt**(1.D0-M))  ! Denominator for the llc calculation
where(abs(deno) > small) llc = llc_ori * (nume / deno)  ! Compute llc where denominator is not too small
where(abs(deno) <= small) llc = llc_ori * 0.5D0         ! Set llc to half of original value if denominator is small

! Calculate alp_s (slope coefficient) and alp_r (river coefficient)
where(qstr_clmt > small) alp_s = (rho**(-M) * qstr_clmt**(M-mm) * Kstr * (0.5D0*lstr)**(-1.D0))**(1.D0/(1.D0-mm))  ! For non-zero streamflow
where(qstr_clmt <= small) alp_s = 0.D0  ! If streamflow is too small, set alp_s to 0

where(qri_clmt > small) alp_r = (rho**(-M) * qri_clmt**(M-mm) * K * llc**(-1.D0))**(1.D0/(1.D0-mm))  ! For non-zero river input
where(qri_clmt <= small) alp_r = 0.D0  ! If river input is too small, set alp_r to 0

!temporal loop
DO iter=step_start,step_end

  ! Read the state of the system from a restart file for the current iteration
  call read_restart(iter,is_coldstart,ny,nc,days_acc_year,days_acc_noleap,days_acc_leap,Ws,Wr,Wr_res,Wr_lake)
  
  ! Read runoff data for the current time step
  call read_runoff(nc,ny,iter,days_acc_year,days_acc_noleap,days_acc_leap,Qrunf,yr_s,mon_s,day_s,d_res,mon_cur)

  !$omp parallel default(shared)
  !$omp workshare

  ! Update state variables: ks, Ws, and Qs 
  where(Qrunf<=small)Qrunf=0.D0  ! Set runoff to zero if it's too small
  Qs0=max(0.D0,alp_s * Ws**(1.D0/(1.D0-mm))) ! Initial flow from stream storage (kg/s)
  ks=max(0.D0,(alp_s/(1.D0-mm)) * Ws**(mm/(1.D0-mm))) ! Flow coefficient (s^-1)
  Ws_last=Ws  ! Store the current water storage 
  where(ks>small) Ws=Ws + (Qrunf-Qs0)/ks*(1.D0-exp(-ks*dt)) ! Update storage (kg)
  where(ks<=small) Ws=Ws + (Qrunf-Qs0)*dt  ! Simplified update if ks is small
  Ws=max(0.D0,Ws)  ! Ensure storage is non-negative
  Qs=max(0.D0,Qrunf-(Ws-Ws_last)/dt)  ! Calculate the stream flow (kg/s)

  ! Calculate variables related to river routing: Qr0, kr
  Qr0=max(0.D0,alp_r * Wr**(1.D0/(1.D0-mm))) ! River flow based on water storage (kg/s)
  kr=max(0.D0,(alp_r/(1.D0-mm)) * Wr**(mm/(1.D0-mm))) ! Flow coefficient for river (s^-1)
  
  ! Update Cl and Al
  where(kr>small.and.abs(kr-ks)>small) Cl=Wr + (Qrunf-Qr0)/kr*(1.D0-exp(-kr*dt)) + (Qrunf-Qs0)/(kr-ks)*(exp(-kr*dt)-exp(-ks*dt))
  where(kr>small.and.abs(kr-ks)<=small) Cl=Wr + (Qrunf-Qr0)/kr*(1.D0-exp(-kr*dt)) - (Qrunf-Qs0)*dt*exp(-kr*dt)
  where(kr<=small.and.ks>small) Cl=Wr + (Qrunf-Qr0)*dt - (Qrunf-Qs0)/ks*(1.D0-exp(-ks*dt))
  where(kr<=small.and.ks<=small) Cl=Wr + (Qs0-Qr0)*dt
  Al=Qs+Wr/dt-Cl/dt  ! Update flow variables

  ! Initialize variables for river routing process
  C1=0.D0
  C2=0.D0
  Qin=0.D0
  Qout=0.D0
  A1=0.D0
  P=0.D0
  B1=0.D0
  B2=0.D0

  !$omp end workshare
  !$omp end parallel

  ! Reservoir module: reset reservoir flow
  Q_res=0.D0
  if(d_res==366)d_res=365  ! Handle leap year day adjustment

  ! Lake module: reset lake flow
  Q_lake=0.D0

  ! Process river routing by going through each node from upstream to downstream
  do n=nmax,0,-1

    !$OMP PARALLEL default(shared) private(i,j,ui,co1,co2,co3)
    !$OMP DO

    ! Loop over each catchment to update the water storage and flow
    do i=1,nc
      if(nts(i)==n)then  ! If the current node matches the iteration step

        ! Process upstream dependencies if any exist
        if(nup(i)>=1)then
          do j=1,nup(i)
            ui=upID(j,i)
            if(ui==-1)exit  ! Exit loop if no more upstream IDs

            ! Calculate flow coefficients based on flow conditions
            if(kr(i)>small)then
              co1=max(0.D0,(1.D0-exp(-kr(i)*dt))/kr(i))
            else
              co1=dt
            endif
            C1(i)=C1(i)+co1*B1(ui)

            if(abs(kr(i)-kr(ui))>small)then
              co2=-(exp(-kr(i)*dt)-exp(-kr(ui)*dt))/(kr(i)-kr(ui))
            else
              co2=dt*exp(-kr(i)*dt)
            endif
            C2(i)=C2(i)+co2*B2(ui)

            ! Process reservoir and lake flows, if active
            if(active_res(ui)==1.and.active_lake(ui)==0)then
              Qin(i)=Qin(i)+Q_res(ui)
            else if(active_res(ui)==0.and.active_lake(ui)==1)then
              Qin(i)=Qin(i)+Q_lake(ui)
            else if(active_res(ui)==1.and.active_lake(ui)==1)then
              Qin(i)=Qin(i)+Q_res(ui)
            else
              Qin(i)=Qin(i)+Qout(ui)
            endif
          enddo
        endif

        ! Update water storage in the current node
        Wr(i)=max(0.D0,Cl(i)+C1(i)+C2(i))
        A1(i)=Qin(i)-C1(i)/dt-C2(i)/dt
        Qout(i)=max(0.D0,Al(i)+A1(i))

        ! Calculate flow parameters based on river flow characteristics
        if(kr(i)>small.and.Qin(i)+Qrunf(i)>small)then
          co3=max(0.D0,(1.D0-exp(-kr(i)*dt))/kr(i))
          P(i)=(dt*Qout(i)-co3*Qr0(i))/((Qin(i)+Qrunf(i))*(dt-co3))
          if(P(i)>0.5D0.and.P(i)<1.5D0)then
            B1(i)=P(i)*(Qin(i)+Qrunf(i))
            B2(i)=-P(i)*(Qin(i)+Qrunf(i))+Qr0(i)
          else
            B1(i)=Qout(i)
            B2(i)=0.D0
          endif
        else
          B1(i)=Qout(i)
          B2(i)=0.D0
          P(i)=-9999.
        endif

        ! Call lake and reservoir calculation subroutines
        call lake_cal(active_lake(i),area_lake(i),Q_lake(i),Wr_lake(i),Qout(i),B1(i),B2(i))
        call res_cal(active_res(i),active_lake(i),Qout(i),Q_lake(i),type_res(i),cat2res(i),&
          Q_res(i),wid_res(i),fld_res(i),Wr_res(i),Qfld_thres(i),cap_res(i),B1(i),B2(i))

      endif
    enddo

    !$OMP END DO
    !$OMP END PARALLEL

  enddo

  ! Write the output for the current time step
  call write_output(nc,yr_s,mon_s,day_s,Qout,Ws,Wr,Q_res,Wr_res,Q_lake,Wr_lake)

ENDDO

end
