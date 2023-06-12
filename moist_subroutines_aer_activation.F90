module moist_subroutines_aer_activation

    use MAPL_ConstantsMod

    implicit none

    public Aer_Activation

    private

    ! Real kind for activation.
    integer,public,parameter :: AER_R4 = SELECTED_REAL_KIND(6,37)
    integer,public,parameter :: AER_R8 = SELECTED_REAL_KIND(15,307)
    integer,public,parameter :: AER_PR = AER_R8

    integer, parameter ::  &
    nsmx_par = 20, npgauss=10 !maximum number of 	
    ! nsmx_par !maximum number of modes allowed

    integer, parameter :: ESMF_MAXSTR = 100

    type, public :: AerProps            
    sequence 
        real, dimension(nsmx_par)  :: num !Num conc m-3
        real, dimension(nsmx_par)  :: dpg !dry Geometric size, m
        real, dimension(nsmx_par)  :: sig  !logarithm (base e) of the dry geometric disp
        real, dimension(nsmx_par)  :: den  !dry density , Kg m-3
        real, dimension(nsmx_par)  :: kap !Hygroscopicity parameter 
        real, dimension(nsmx_par)  :: fdust! mass fraction of dust 
        real, dimension(nsmx_par)  :: fsoot ! mass fraction of soot
        real, dimension(nsmx_par)  :: forg ! mass fraction of organics
        integer   :: nmods  ! total number of modes (nmods<nmodmax)
    end type AerProps 

    LOGICAL  :: USE_AEROSOL_NN = .TRUE.

    real        , parameter :: R_AIR     =  3.47e-3 !m3 Pa kg-1K-1
    real(AER_PR), parameter :: zero_par  =  1.e-6   ! small non-zero value
    real(AER_PR), parameter :: ai        =  0.0000594
    real(AER_PR), parameter :: bi        =  3.33
    real(AER_PR), parameter :: ci        =  0.0264
    real(AER_PR), parameter :: di        =  0.0033

    real(AER_PR), parameter :: betai     = -2.262e+3
    real(AER_PR), parameter :: gamai     =  5.113e+6
    real(AER_PR), parameter :: deltai    =  2.809e+3
    real(AER_PR), parameter :: densic    =  917.0   !Ice crystal density in kgm-3

    real, parameter :: NN_MIN      =  100.0e6
    real, parameter :: NN_MAX      = 1000.0e6
    real :: tStart, tEnd, time

    contains

    SUBROUTINE Aer_Activation(IM,JM,LM, q, t, plo, ple, zlo, zle, qlcn, qicn, qlls, qils, &
            sh, evap, kpbl, tke, vvel, FRLAND, USE_AERO_BUFFER, &
            AeroProps, aero_aci, NACTL, NACTI, NWFA, NN_LAND, NN_OCEAN,dirName, rank_str)
        IMPLICIT NONE
        integer, intent(in)::IM,JM,LM
        TYPE(AerProps), dimension (IM,JM,LM),intent(inout)  :: AeroProps
        ! type(ESMF_State)            ,intent(inout) :: aero_aci
        integer, intent(inout) :: aero_aci
        real, dimension (IM,JM,LM)  ,intent(in ) :: plo ! Pa
        real, dimension (IM,JM,0:LM),intent(in ) :: ple ! Pa
        real, dimension (IM,JM,LM)  ,intent(in ) :: q,t,tke,vvel,zlo, qlcn, qicn, qlls, qils
        real, dimension (IM,JM,0:LM),intent(in ) :: zle
        real, dimension (IM,JM)     ,intent(in ) :: FRLAND
        real, dimension (IM,JM)     ,intent(in ) :: sh, evap, kpbl
        real                        ,intent(in ) :: NN_LAND, NN_OCEAN 
        logical                     ,intent(in ) :: USE_AERO_BUFFER

        real, dimension (IM,JM,LM),intent(OUT) :: NACTL,NACTI, NWFA

        real(AER_PR), allocatable, dimension (:) :: sig0,rg,ni,bibar,nact 
        real(AER_PR)                     :: wupdraft,tk,press,air_den,QI,QL,WC,BB,RAUX

        integer, dimension (IM,JM) :: kpbli

        real, dimension(:,:,:,:,:), allocatable :: buffer

        character(len=ESMF_MAXSTR)      :: aci_field_name
        real, pointer, dimension(:,:)   :: aci_ptr_2d
        real, pointer, dimension(:,:,:) :: aci_ptr_3d
        real, pointer, dimension(:,:,:) :: aci_num
        real, pointer, dimension(:,:,:) :: aci_dgn
        real, pointer, dimension(:,:,:) :: aci_sigma
        real, pointer, dimension(:,:,:) :: aci_density
        real, pointer, dimension(:,:,:) :: aci_hygroscopicity
        real, pointer, dimension(:,:,:) :: aci_f_dust
        real, pointer, dimension(:,:,:) :: aci_f_soot
        real, pointer, dimension(:,:,:) :: aci_f_organic
        character(len=ESMF_MAXSTR), allocatable, dimension(:) :: aero_aci_modes
        integer                         :: ACI_STATUS
        REAL                            :: aux1,aux2,aux3,hfs,hfl,nfaux

        integer :: n_modes
        REAL :: numbinit
        integer :: i,j,k,n,rc

        character(len=ESMF_MAXSTR)              :: IAm="Aer_Activation"
        integer                                 :: STATUS

        character(2) :: n_str
        character*100 :: dirName, rank_str
        integer :: rank, fileID

        do k = 1, LM
        do j = 1, JM
        do i = 1, IM
        AeroProps(i,j,k)%num = 0.0
        end do
        end do
        end do    
        NACTL    = 0.
        NACTI    = 0.

        kpbli = MAX(MIN(NINT(kpbl),LM-1),1)

        if (USE_AEROSOL_NN) then

        ! call ESMF_AttributeGet(aero_aci, name='number_of_aerosol_modes', value=n_modes, __RC__)

            open(newunit=fileID, file=trim(dirName) // "/n_modes_" // trim(rank_str) // ".in", &
                status='old', form="unformatted", action="read")
            read(fileID) n_modes
            close(fileID)

            if (n_modes > 0) then

        ! allocate( sig0(n_modes), __STAT__)
        ! allocate(   rg(n_modes), __STAT__)
        ! allocate(   ni(n_modes), __STAT__)
        ! allocate(bibar(n_modes), __STAT__)
        ! allocate( nact(n_modes), __STAT__)

        ! allocate(aero_aci_modes(n_modes), __STAT__)
        ! call ESMF_AttributeGet(aero_aci, name='aerosol_modes', itemcount=n_modes, valuelist=aero_aci_modes, __RC__)

        ! call ESMF_AttributeGet(aero_aci, name='air_pressure_for_aerosol_optics', value=aci_field_name, __RC__)
        ! if (aci_field_name /= '') then
        ! call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
        ! aci_ptr_3d = PLE
        ! end if

        ! call ESMF_AttributeGet(aero_aci, name='air_temperature', value=aci_field_name, __RC__)
        ! if (aci_field_name /= '') then
        ! call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
        ! aci_ptr_3d = T
        ! end if

        ! call ESMF_AttributeGet(aero_aci, name='fraction_of_land_type', value=aci_field_name, __RC__)
        ! if (aci_field_name /= '') then
        ! call MAPL_GetPointer(aero_aci, aci_ptr_2d, trim(aci_field_name), __RC__)
        ! aci_ptr_2d = FRLAND
        ! end if

                if (USE_AERO_BUFFER) then
                ! allocate(buffer(im,jm,lm,n_modes,8), __STAT__)
                    allocate(buffer(im,jm,lm,n_modes,8))
                end if

                allocate( sig0(n_modes))
                allocate(   rg(n_modes))
                allocate(   ni(n_modes))
                allocate(bibar(n_modes))
                allocate( nact(n_modes))
                ! allocate(ac(n_modes))
                ! allocate(fracactn(n_modes))
                ! allocate(eta(n_modes))
                ! allocate(xlogsigm(n_modes))
                ! allocate(sm(n_modes))

                allocate(aero_aci_modes(n_modes))

                ! *** Allocations done for standalone ***
                allocate(aci_num            (im, jm, lm))
                allocate(aci_dgn            (im, jm, lm))
                allocate(aci_sigma          (im, jm, lm))
                allocate(aci_density        (im, jm, lm))
                allocate(aci_hygroscopicity (im, jm, lm))
                allocate(aci_f_dust         (im, jm, lm))
                allocate(aci_f_soot         (im, jm, lm))
                allocate(aci_f_organic      (im, jm, lm))

                ACTIVATION_PROPERTIES: do n = 1, n_modes


                    if(n.ge.10) then
                        write(n_str,'(i2)') n
                    else
                        write(n_str,'(i1)') n
                    endif

        ! call ESMF_AttributeSet(aero_aci, name='aerosol_mode', value=trim(aero_aci_modes(n)), __RC__)

        ! ! execute the aerosol activation properties method 
        ! call ESMF_MethodExecute(aero_aci, label='aerosol_activation_properties', userRC=ACI_STATUS, RC=STATUS)
        ! VERIFY_(ACI_STATUS)
        ! VERIFY_(STATUS)

        ! ! copy out aerosol activation properties
        ! call ESMF_AttributeGet(aero_aci, name='aerosol_number_concentration', value=aci_field_name, __RC__)
        ! call MAPL_GetPointer(aero_aci, aci_num, trim(aci_field_name), __RC__)

        ! call ESMF_AttributeGet(aero_aci, name='aerosol_dry_size', value=aci_field_name, __RC__)
        ! call MAPL_GetPointer(aero_aci, aci_dgn, trim(aci_field_name), __RC__)

        ! call ESMF_AttributeGet(aero_aci, name='width_of_aerosol_mode', value=aci_field_name, __RC__)
        ! call MAPL_GetPointer(aero_aci, aci_sigma, trim(aci_field_name), __RC__)

        ! call ESMF_AttributeGet(aero_aci, name='aerosol_density', value=aci_field_name, __RC__)
        ! call MAPL_GetPointer(aero_aci, aci_density, trim(aci_field_name), __RC__)

        ! call ESMF_AttributeGet(aero_aci, name='aerosol_hygroscopicity', value=aci_field_name, __RC__)
        ! call MAPL_GetPointer(aero_aci, aci_hygroscopicity, trim(aci_field_name), __RC__)

        ! call ESMF_AttributeGet(aero_aci, name='fraction_of_dust_aerosol', value=aci_field_name, __RC__)
        ! call MAPL_GetPointer(aero_aci, aci_f_dust, trim(aci_field_name), __RC__)

        ! call ESMF_AttributeGet(aero_aci, name='fraction_of_soot_aerosol', value=aci_field_name, __RC__)
        ! call MAPL_GetPointer(aero_aci, aci_f_soot, trim(aci_field_name), __RC__)

        ! call ESMF_AttributeGet(aero_aci, name='fraction_of_organic_aerosol', value=aci_field_name, __RC__)
        ! call MAPL_GetPointer(aero_aci, aci_f_organic, trim(aci_field_name), __RC__)

                    open(newunit=fileID, file=trim(dirName) // "/aci_num_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
                        status='old', form="unformatted", action="read")
                    read(fileID) aci_num
                    close(fileID)

                    open(newunit=fileID, file=trim(dirName) // "/aci_dgn_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
                        status='old', form="unformatted", action="read")
                    read(fileID) aci_dgn
                    close(fileID)

                    open(newunit=fileID, file=trim(dirName) // "/aci_sigma_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
                        status='old', form="unformatted", action="read")
                    read(fileID) aci_sigma
                    close(fileID)

                    open(newunit=fileID, file=trim(dirName) // "/aci_hygroscopicity_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
                        status='old', form="unformatted", action="read")
                    read(fileID) aci_hygroscopicity
                    close(fileID)

                    open(newunit=fileID, file=trim(dirName) // "/aci_density_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
                        status='old', form="unformatted", action="read")
                    read(fileID) aci_density
                    close(fileID)

                    open(newunit=fileID, file=trim(dirName) // "/aci_f_dust_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
                        status='old', form="unformatted", action="read")
                    read(fileID) aci_f_dust
                    close(fileID)

                    open(newunit=fileID, file=trim(dirName) // "/aci_f_soot_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
                        status='old', form="unformatted", action="read")
                    read(fileID) aci_f_soot
                    close(fileID)

                    open(newunit=fileID, file=trim(dirName) // "/aci_f_organic_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
                        status='old', form="unformatted", action="read")
                    read(fileID) aci_f_organic
                    close(fileID)

                    if (USE_AERO_BUFFER) then
                        buffer(:,:,:,n,1) = aci_num
                        buffer(:,:,:,n,2) = aci_dgn
                        buffer(:,:,:,n,3) = aci_sigma
                        buffer(:,:,:,n,4) = aci_hygroscopicity
                        buffer(:,:,:,n,5) = aci_density
                        buffer(:,:,:,n,6) = aci_f_dust
                        buffer(:,:,:,n,7) = aci_f_soot
                        buffer(:,:,:,n,8) = aci_f_organic
                    else
                        AeroProps(:,:,:)%num(n)   = aci_num
                        AeroProps(:,:,:)%dpg(n)   = aci_dgn
                        AeroProps(:,:,:)%sig(n)   = aci_sigma
                        AeroProps(:,:,:)%kap(n)   = aci_hygroscopicity
                        AeroProps(:,:,:)%den(n)   = aci_density
                        AeroProps(:,:,:)%fdust(n) = aci_f_dust
                        AeroProps(:,:,:)%fsoot(n) = aci_f_soot
                        AeroProps(:,:,:)%forg(n)  = aci_f_organic
                        AeroProps(:,:,:)%nmods    = n_modes                 ! no need of a 3D field: aero provider specific
                    end if

                end do ACTIVATION_PROPERTIES

                call cpu_time(tStart)

!$acc data present(AeroProps,AeroProps%num,AeroProps%dpg,AeroProps%sig,AeroProps%kap,AeroProps%den,AeroProps%fdust,&
!$acc              AeroProps%fsoot,AeroProps%forg,&
!$acc              NWFA,ple,t,q,sh,evap,zle,plo,qicn,qils,qlcn,qlls,tke, vvel,NACTL,NACTI,FRLAND) &
!$acc      copyin(kpbli,buffer) create(ni,rg,sig0,nact)

                if (USE_AERO_BUFFER) then
!$acc parallel loop gang vector collapse(4)
                    do k = 1, LM
                        do j = 1, JM
                            do i = 1, IM
                                do n = 1, n_modes
                                    AeroProps(i,j,k)%num(n)   = buffer(i,j,k,n,1)
                                    AeroProps(i,j,k)%dpg(n)   = buffer(i,j,k,n,2)
                                    AeroProps(i,j,k)%sig(n)   = buffer(i,j,k,n,3)
                                    AeroProps(i,j,k)%kap(n)   = buffer(i,j,k,n,4)
                                    AeroProps(i,j,k)%den(n)   = buffer(i,j,k,n,5)
                                    AeroProps(i,j,k)%fdust(n) = buffer(i,j,k,n,6)
                                    AeroProps(i,j,k)%fsoot(n) = buffer(i,j,k,n,7)
                                    AeroProps(i,j,k)%forg(n)  = buffer(i,j,k,n,8)
                                end do
                                AeroProps(i,j,k)%nmods       = n_modes                 ! no need of a 3D field: aero provider specific
                            end do
                        end do
                    end do
!$acc end parallel
                    deallocate(buffer)
                end if
                 
!$acc parallel loop gang vector collapse(3) private(nfaux)
                do k = 1, LM
                    do j = 1, JM
                        do i = 1, IM
                        nfaux =  0.0
!$acc loop seq
                        do n = 1, n_modes
                            if (AeroProps(i,j,k)%kap(n) .gt. 0.4) then 
                                    nfaux =  nfaux + AeroProps(i,j,k)%num(n)
                            end if                           
                        end do !modes
                        NWFA(I, J, K)  =  nfaux
                        end do
                    end do 
                end do 
!$acc end parallel
   
                deallocate(aero_aci_modes)
   
                !--- activated aerosol # concentration for liq/ice phases (units: m^-3)
                numbinit = 0.
                WC       = 0.
                BB       = 0.
                RAUX     = 0.

!$acc parallel loop gang vector collapse(3) &
!$acc          private(tk, press, air_den, qi, ql, wupdraft, numbinit)
                DO k=LM,1,-1
                    DO j=1,JM
                        DO i=1,IM
                            NACTL(i,j,k) = NN_LAND*FRLAND(i,j) + NN_OCEAN*(1.0-FRLAND(i,j))
                            NACTI(i,j,k) = NN_LAND*FRLAND(i,j) + NN_OCEAN*(1.0-FRLAND(i,j))
                        
                            tk                 = T(i,j,k)                         ! K
                            press              = plo(i,j,k)                       ! Pa   
                            air_den            = press*28.8e-3/8.31/tk            ! kg/m3
                            qi                 = (qicn(i,j,k)+qils(i,j,k))*1.e+3  ! g/kg
                            ql                 = (qlcn(i,j,k)+qlls(i,j,k))*1.e+3  ! g/kg
                            wupdraft           = vvel(i,j,k) + SQRT(tke(i,j,k))
            
                            ! Liquid Clouds
                            IF( (tk >= MAPL_TICE-40.0) .and. (plo(i,j,k) > 10000.0) .and. &
                                (wupdraft > 0.1 .and. wupdraft < 100.) ) then
                
                                ni   (1:n_modes)    =   max(AeroProps(i,j,k)%num(1:n_modes)*air_den,  zero_par)  ! unit: [m-3]
                                rg   (1:n_modes)    =   max(AeroProps(i,j,k)%dpg(1:n_modes)*0.5*1.e6, zero_par)  ! unit: [um]
                                sig0 (1:n_modes)    =       AeroProps(i,j,k)%sig(1:n_modes)                      ! unit: [um]
                                bibar(1:n_modes)    =   max(AeroProps(i,j,k)%kap(1:n_modes),          zero_par)                 
                
                                call GetActFrac(           n_modes    &
                                                ,      ni(1:n_modes)   &
                                                ,      rg(1:n_modes)   &
                                                ,    sig0(1:n_modes)   &
                                                ,      tk              &
                                                ,   press              &
                                                ,wupdraft              &
                                                ,    nact(1:n_modes)   &
                                                ,   bibar(1:n_modes)   &
                                                )
        
                                numbinit = 0.
                                NACTL(i,j,k) = 0.
                                DO n=1,n_modes
                                    numbinit = numbinit + AeroProps(i,j,k)%num(n)*air_den
                                    NACTL(i,j,k)= NACTL(i,j,k) + nact(n) !#/m3
                                ENDDO
                                NACTL(i,j,k) = MIN(NACTL(i,j,k),0.99*numbinit)
        
                            ENDIF ! Liquid Clouds
        
                            ! Ice Clouds
                            IF( (tk <= MAPL_TICE) .and. ((QI > tiny(1.)) .or. (QL > tiny(1.))) ) then
                                numbinit = 0.
!$acc loop seq
                                DO n=1,n_modes
                                    if (AeroProps(i,j,k)%dpg(n) .ge. 0.5e-6) & ! diameters > 0.5 microns
                                    numbinit = numbinit + AeroProps(i,j,k)%num(n)
                                ENDDO
                                numbinit = numbinit * air_den ! #/m3
                                ! Number of activated IN following deMott (2010) [#/m3]
                                NACTI(i,j,k) = ai*((MAPL_TICE-tk)**bi) * numbinit**(ci*(MAPL_TICE-tk)+di) !#/m3
                            ENDIF
        
                            !-- apply limits for NACTL/NACTI
                            IF(NACTL(i,j,k) < NN_MIN) NACTL(i,j,k) = NN_MIN
                            IF(NACTL(i,j,k) > NN_MAX) NACTL(i,j,k) = NN_MAX
                            IF(NACTI(i,j,k) < NN_MIN) NACTI(i,j,k) = NN_MIN
                            IF(NACTI(i,j,k) > NN_MAX) NACTI(i,j,k) = NN_MAX
        
                ENDDO;ENDDO;ENDDO
!$acc end data 
            deallocate(   rg)
            deallocate(   ni)
            deallocate(bibar)
            deallocate( nact)

        end if ! n_modes > 0
         
        else ! USE_AEROSOL_NN
!$acc kernels
            do k = 1, LM
                NACTL(:,:,k) = NN_LAND*FRLAND + NN_OCEAN*(1.0-FRLAND)
                NACTI(:,:,k) = NN_LAND*FRLAND + NN_OCEAN*(1.0-FRLAND)
            end do
!$acc end kernels
        end if



                call cpu_time(tEnd)
                time = tEnd - tStart
                print*, 'Aer_Activation Computation Time = ', time

    END SUBROUTINE Aer_Activation

    subroutine GetActFrac(nmodes  & !nmodes                                     &
        ,xnap                         & !ni              (1:nmodes)              & 
        ,rg                           & !0.5d+00*dgn_dry (1:nmodes)              & 
        ,sigmag                       & !sig0             (1:nmodes)              & 
        ,tkelvin                      & !tk              (i,j,k)                     &
        ,ptot                         & !pres             (i,j,k)                     & 
        ,wupdraft                     & !wupdraft             (i,j,k)                     & 
        ,nact                         & !nact             (i,j,k,1:nmodes)             &
        ,bibar)
!$acc routine seq    
        IMPLICIT NONE
  
        ! arguments.
        
        integer :: nmodes               !< number of modes [1]      
        real(AER_PR) :: xnap(nmodes)         !< number concentration for each mode [#/m^3]
        real(AER_PR) :: rg(nmodes)           !< geometric mean dry radius for each mode [um]
        real(AER_PR) :: sigmag(nmodes)       !< geometric standard deviation for each mode [um]
        real(AER_PR) :: tkelvin              !< absolute temperature [k]
        real(AER_PR) :: ptot                 !< ambient pressure [pa]
        real(AER_PR) :: wupdraft             !< updraft velocity [m/s]
!       real(AER_PR) :: ac(nmodes)           !< minimum dry radius for activation for each mode [um]
!       real(AER_PR) :: fracactn(nmodes)     !< activating fraction of number conc. for each mode [1]
        real(AER_PR) :: nact(nmodes)         !< activating number concentration for each mode [#/m^3]
        real(AER_PR) :: bibar(nmodes)        ! hygroscopicity parameter for each mode [1]
  
        ! local variables. 
  
        integer :: i, j                 ! loop counters       
  
        !--------------------------------------------------------------------------------------------------------------
        ! calculate the droplet activation parameters for each mode. 
        !--------------------------------------------------------------------------------------------------------------
        call ActFrac_Mat(nmodes,xnap,rg,sigmag,bibar,tkelvin,ptot,wupdraft,nact)
  
    end subroutine GetActFrac

    subroutine ActFrac_Mat(nmodes,xnap,rg,sigmag,bibar,tkelvin,ptot,wupdraft,nact)
!$acc routine seq
        IMPLICIT NONE
  
        ! Arguments.
        
        integer :: nmodes            !< number of modes [1]      
        real(AER_PR) :: xnap(nmodes)      !< number concentration for each mode [#/m^3]
  !     real(AER_PR) :: xmap(nmodes)      !< mass   concentration for each mode [ug/m^3]
        real(AER_PR) :: rg(nmodes)        !< geometric mean radius for each mode [um]
        real(AER_PR) :: sigmag(nmodes)    !< geometric standard deviation for each mode [um]
        real(AER_PR) :: bibar(nmodes)     !< hygroscopicity parameter for each mode [1]
        real(AER_PR) :: tkelvin           !< absolute temperature [k]
        real(AER_PR) :: ptot              !< ambient pressure [pa]
        real(AER_PR) :: wupdraft          !< updraft velocity [m/s]
        ! real(AER_PR) :: ac(nmodes)        !< minimum dry radius for activation for each mode [um]
        ! real(AER_PR) :: fracactn(nmodes)  !< activating fraction of number conc. for each mode [1]
        real(AER_PR) :: nact(nmodes)      !< activating number concentration for each mode [#/m^3]
  
        ! parameters.
        
        real(AER_PR), parameter :: pi            = 3.141592653589793d+00
        real(AER_PR), parameter :: twopi         = 2.0d+00 * pi
        real(AER_PR), parameter :: sqrt2         = 1.414213562d+00
        real(AER_PR), parameter :: threesqrt2by2 = 1.5d+00 * sqrt2
  
        real(AER_PR), parameter :: avgnum   = 6.0221367d+23       ! [1/mol]
        real(AER_PR), parameter :: rgasjmol = 8.31451d+00         ! [j/mol/k]
        real(AER_PR), parameter :: wmolmass = 18.01528d-03        ! molar mass of h2o     [kg/mol]
        real(AER_PR), parameter :: amolmass = 28.966d-03          ! molar mass of air     [kg/mol]
        real(AER_PR), parameter :: asmolmss = 132.1406d-03        ! molar mass of nh42so4 [kg/mol]
        real(AER_PR), parameter :: denh2o   = 1.00d+03            ! density of water [kg/m^3]
        real(AER_PR), parameter :: denamsul = 1.77d+03            ! density of pure ammonium sulfate [kg/m^3]
        real(AER_PR), parameter :: xnuamsul = 3.00d+00            ! # of ions formed when the salt is dissolved in water [1]
        real(AER_PR), parameter :: phiamsul = 1.000d+00           ! osmotic coefficient value in a-r 1998. [1] 
        real(AER_PR), parameter :: gravity  = 9.81d+00            ! grav. accel. at the earth's surface [m/s/s] 
        real(AER_PR), parameter :: heatvap  = 40.66d+03/wmolmass  ! latent heat of vap. for water and tnbp [j/kg] 
        real(AER_PR), parameter :: cpair    = 1006.0d+00          ! heat capacity of air [j/kg/k] 
        real(AER_PR), parameter :: t0dij    = 273.15d+00          ! reference temp. for dv [k] 
        real(AER_PR), parameter :: p0dij    = 101325.0d+00        ! reference pressure for dv [pa] 
        real(AER_PR), parameter :: dijh2o0  = 0.211d-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
        !----------------------------------------------------------------------------------------------------------------    
        ! real(AER_PR), parameter :: t0dij    = 283.15d+00          ! reference temp. for dv [k] 
        ! real(AER_PR), parameter :: p0dij    = 80000.0d+00         ! reference pressure for dv [pa] 
        ! real(AER_PR), parameter :: dijh2o0  = 0.300d-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
        !----------------------------------------------------------------------------------------------------------------
        real(AER_PR), parameter :: deltav   = 1.096d-07           ! vapor jump length [m]  
        real(AER_PR), parameter :: deltat   = 2.160d-07           ! thermal jump length [m]  
        real(AER_PR), parameter :: alphac   = 1.000d+00           ! condensation mass accommodation coefficient [1]  
        real(AER_PR), parameter :: alphat   = 0.960d+00           ! thermal accommodation coefficient [1]  
  
        ! local variables. 
  
        integer            :: i                              ! loop counter 
        real(AER_PR)            :: dv                             ! diffusion coefficient for water [m^2/s] 
        real(AER_PR)            :: dvprime                        ! modified diffusion coefficient for water [m^2/s] 
        real(AER_PR)            :: dumw, duma                     ! scratch variables [s/m] 
        real(AER_PR)            :: wpe                            ! saturation vapor pressure of water [pa]  
        real(AER_PR)            :: surten                         ! surface tension of air-water interface [j/m^2] 
        real(AER_PR)            :: xka                            ! thermal conductivity of air [j/m/s/k]  
        real(AER_PR)            :: xkaprime                       ! modified thermal conductivity of air [j/m/s/k]  
        real(AER_PR)            :: eta                            ! model parameter [1]  
        real(AER_PR)            :: zeta                           ! model parameter [1]  
        ! real(AER_PR)            :: xlogsigm(nmodes)               ! ln(sigmag) [1]   
        real(AER_PR)            :: a                              ! [m]
        real(AER_PR)            :: g                              ! [m^2/s]   
        real(AER_PR)            :: rdrp                           ! [m]   
        real(AER_PR)            :: f1                             ! [1]   
        real(AER_PR)            :: f2                             ! [1]
        real(AER_PR)            :: alpha                          ! [1/m]
        real(AER_PR)            :: gamma                          ! [m^3/kg]   
        real(AER_PR)            :: sm                             ! [1]   
        real(AER_PR)            :: dum                            ! [1/m]    
        real(AER_PR)            :: u                              ! argument to error function [1]
        real(AER_PR)            :: erf                            ! error function [1], but not declared in an f90 module 
        real(AER_PR)            :: smax                           ! maximum supersaturation [1]
  
        real(AER_PR)            :: aux1,aux2                      ! aux
        real(AER_PR)            :: ac                             !< minimum dry radius for activation for each mode [um]
        real(AER_PR)            :: fracactn                       !< activating fraction of number conc. for each mode [1]
  !----------------------------------------------------------------------------------------------------------------------
  !     rdrp is the radius value used in eqs.(17) & (18) and was adjusted to yield eta and zeta 
  !     values close to those given in a-z et al. 1998 figure 5. 
  !----------------------------------------------------------------------------------------------------------------------
        rdrp = 0.105d-06   ! [m] tuned to approximate the results in figures 1-5 in a-z et al. 1998.  
  !----------------------------------------------------------------------------------------------------------------------
  !     these variables are common to all modes and need only be computed once. 
  !----------------------------------------------------------------------------------------------------------------------
        dv = dijh2o0*(p0dij/ptot)*(tkelvin/t0dij)**1.94d+00                 ! [m^2/s] (p&k,2nd ed., p.503)
        surten = 76.10d-03 - 0.155d-03 * (tkelvin-273.15d+00)               ! [j/m^2] 
        wpe = exp( 77.34491296d+00 - 7235.424651d+00/tkelvin - 8.2d+00*log(tkelvin) + tkelvin*5.7113d-03 )  ! [pa] 
        dumw = sqrt(twopi*wmolmass/rgasjmol/tkelvin)                        ! [s/m] 
        dvprime = dv / ( (rdrp/(rdrp+deltav)) + (dv*dumw/(rdrp*alphac)) )   ! [m^2/s] - eq. (17) 
        xka = (5.69d+00+0.017d+00*(tkelvin-273.15d+00))*418.4d-05           ! [j/m/s/k] (0.0238 j/m/s/k at 273.15 k)
        duma = sqrt(twopi*amolmass/rgasjmol/tkelvin)                        ! [s/m]
        xkaprime = xka / ( ( rdrp/(rdrp+deltat) ) + ( xka*duma/(rdrp*alphat*denh2o*cpair) ) )   ! [j/m/s/k]
        g = 1.0d+00 / ( (denh2o*rgasjmol*tkelvin) / (wpe*dvprime*wmolmass) &
                        + ( (heatvap*denh2o) / (xkaprime*tkelvin) ) &
                        * ( (heatvap*wmolmass) / (rgasjmol*tkelvin) - 1.0d+00 ) )               ! [m^2/s]
        a = (2.0d+00*surten*wmolmass)/(denh2o*rgasjmol*tkelvin)                                 ! [m] 
        alpha = (gravity/(rgasjmol*tkelvin))*((wmolmass*heatvap)/(cpair*tkelvin) - amolmass)    ! [1/m] 
        gamma = (rgasjmol*tkelvin)/(wpe*wmolmass) &
              + (wmolmass*heatvap*heatvap)/(cpair*ptot*amolmass*tkelvin)                        ! [m^3/kg]
        dum = sqrt(alpha*wupdraft/g)                  ! [1/m] 
        zeta = 2.d+00*a*dum/3.d+00                    ! [1] 
  !----------------------------------------------------------------------------------------------------------------
        ! write(1,'(a27,4d15.5)')'surten,wpe,a            =',surten,wpe,a
        ! write(1,'(a27,4d15.5)')'xka,xkaprime,dv,dvprime =',xka,xkaprime,dv,dvprime
        ! write(1,'(a27,4d15.5)')'alpha,gamma,g, zeta     =',alpha,gamma,g,zeta
  !----------------------------------------------------------------------------------------------------------------------
  !     these variables must be computed for each mode. 
  !----------------------------------------------------------------------------------------------------------------------
        ! xlogsigm(:) = log(sigmag(:))                                                   ! [1] 
        smax = 0.0d+00                                                                   ! [1]
        
        do i=1, nmodes
      
            sm    = ( 2.0d+00/sqrt(bibar(i)) ) * ( a/(3.0*rg(i)) )**1.5d+00              ! [1] 
            eta   = dum**3 / (twopi*denh2o*gamma*xnap(i))                                ! [1] 
    
            !--------------------------------------------------------------------------------------------------------------
            ! write(1,'(a27,i4,4d15.5)')'i,eta(i),sm(i) =',i,eta(i),sm(i)
            !--------------------------------------------------------------------------------------------------------------
            f1 = 0.5d+00 * exp(2.50d+00 * log(sigmag(i))**2)                              ! [1] 
            f2 = 1.0d+00 +     0.25d+00 * log(sigmag(i))                                  ! [1] 
            smax = smax + (   f1*(  zeta  / eta                 )**1.50d+00 &
                            + f2*(sm**2/(eta   +3.0d+00*zeta))**0.75d+00 ) / sm**2        ! [1] - eq. (6)
        enddo 
        smax = 1.0d+00 / sqrt(smax)                                                       ! [1]
        
        
        !aux1=0.0D0; aux2=0.0D0
        do i=1, nmodes
            sm    = ( 2.0d+00/sqrt(bibar(i)) ) * ( a/(3.0*rg(i)) )**1.5d+00               ! [1] 
            ac          = rg(i) * ( sm / smax )**0.66666666666666667d+00               ! [um]
    
            u           = log(ac/rg(i)) / ( sqrt2 * log(sigmag(i)) )                      ! [1]
            fracactn    = 0.5d+00 * (1.0d+00 - erf(u))                                    ! [1]
            nact(i)     = fracactn * xnap(i)                                              ! [#/m^3]
          !--------------------------------------------------------------------------------------------------------------
          !aux1=aux1+ xnap(i); aux2=aux2+ nact(i) 
          
          !if(fracactn(i) .gt. 0.9999999d+00 ) then
          !   write(*,*)i,ac(i),u,fracactn(i),xnap(i)
          !  print*,' xxx',i,ac(i),u,fracactn(i),xnap(i)
          ! stop
          ! endif
          !--------------------------------------------------------------------------------------------------------------
        enddo 
  
        return
    end subroutine ActFrac_Mat
end module

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”
  
! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”
  
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at
  
! http://www.apache.org/licenses/LICENSE-2.0
  
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.
