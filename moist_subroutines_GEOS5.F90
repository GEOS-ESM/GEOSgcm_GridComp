module moist_subroutines_GEOS5

    use MAPL_ConstantsMod
    use ConvPar_GF_SharedParams
    use module_gate

    implicit none

    private

    public CUP_gf_sh, CUP_gf_GEOS5

    ! Real kind for activation.
    integer,public,parameter :: AER_R4 = SELECTED_REAL_KIND(6,37)
    integer,public,parameter :: AER_R8 = SELECTED_REAL_KIND(15,307)
    integer,public,parameter :: AER_PR = AER_R8

    real        , parameter :: R_AIR     =  3.47e-3 !m3 Pa kg-1K-1

    ! combined constantc
    real, parameter :: cpbgrav = MAPL_CP/MAPL_GRAV
    real, parameter :: gravbcp = MAPL_GRAV/MAPL_CP
    real, parameter :: alhlbcp = MAPL_ALHL/MAPL_CP
    real, parameter :: alhfbcp = MAPL_ALHF/MAPL_CP
    real, parameter :: alhsbcp = MAPL_ALHS/MAPL_CP

    CHARACTER(LEN=10),PARAMETER  :: host_model = 'NEW_GEOS5'

    !-autonversion formulation: (1) Kessler, (2) Berry, (3) NOAA, (4) Sundvisqt
    INTEGER, PARAMETER :: autoconv = 1

    !- proportionality constant to estimate pressure
    !- gradient of updraft (Zhang and Wu, 2003, JAS)
    REAL, PARAMETER ::    pgcd=1., pgcon=0.

    ! (do not change  to false)

    !- numerical constraints
    REAL, PARAMETER ::       &
        xmbmaxshal  =  0.05,    &  ! kg/m2/s
        mintracer   =  tiny(1.),&  ! kg/kg - tiny(x)
        smallerQV   =  1.e-16      ! kg/kg  

    !- for internal debugging
    LOGICAL :: wrtgrads = .false.

 !-- General control for the diverse options in GF
    LOGICAL, PARAMETER :: ENTRNEW        = .TRUE.  !- new entr formulation
    LOGICAL, PARAMETER :: COUPL_MPHYSICS = .TRUE.  !- coupling with cloud microphysics
                                                   ! (do not change  to false)
    LOGICAL, PARAMETER :: MELT_GLAC      = .TRUE.  !- turn ON/OFF ice phase/melting
    LOGICAL, PARAMETER :: DOWNDRAFT      = .TRUE.  !- turn ON/OFF downdrafts
    LOGICAL, PARAMETER :: MOMENTUM       = .TRUE.  !- turn ON/OFF conv transp of momentum
   
    LOGICAL, PARAMETER :: FEED_3DMODEL   = .TRUE.  !- set "false" to not feedback the AGCM with the
                                                   !- heating/drying/transport conv tendencies
    LOGICAL, PARAMETER :: USE_C1D        = .TRUE.  !- turn ON/OFF the 'c1d' detrainment approach
   
    INTEGER, PARAMETER :: VERT_DISCR = 1
   
    !-autonversion formulation: (1) original , (2) MP_GT
    INTEGER, PARAMETER :: CLOUDMP = 1
   
    !-rainfall evaporation(1) orig (2) mix orig+new (3) new
    INTEGER, PARAMETER :: aeroevap = 1
   
    INTEGER, PARAMETER ::               &
                          maxens  = 1,  & ! 1  ensemble one on cap_max
                          maxens2 = 1,  & ! 1  ensemble two on precip efficiency
                          maxens3 = 16, & !16 ensemble three done in cup_forcing_ens16 for G3d
                          ensdim  = maxens*maxens2*maxens3,&
                          ens4    = 1

    ! !- physical constants
     REAL, PARAMETER ::  &
        ! rgas    = 287.,    & ! J K-1 kg-1
        ! cp      = 1004.,   & ! J K-1 kg-1
        ! rv      = 461.,    & ! J K-1 kg-1
        ! p00     = 1.e5,    & ! hPa
        ! tcrit   = 258.,    & ! K
        ! g       = MAPL_GRAV,& ! m s-2
        ! cpor    = cp/rgas, &
        ! xlv     = 2.5e6,   & ! J kg-1
        ! akmin   = 1.0,     & ! #
        ! tkmin   = 1.e-5,   & ! m+2 s-2
        ! ccnclean= 250.,    & ! # cm-3
        ! T_0     = 273.16,  & ! K
        T_ice   = 250.16 !,  & ! K
        ! xlf     = 0.333e6    ! latent heat of freezing (J K-1 kg-1)
    contains

    ! SUBROUTINE Aer_Activation(IM,JM,LM, q, t, plo, ple, zlo, zle, qlcn, qicn, qlls, qils, &
    !     sh, evap, kpbl, omega, FRLAND, USE_AERO_BUFFER, &
    !     AeroProps, aero_aci, NACTL, NACTI, NWFA, dirName, rank_str)
    !     IMPLICIT NONE
    !     integer, intent(in)::IM,JM,LM
    !     TYPE(AerProps), dimension (IM,JM,LM),intent(inout)  :: AeroProps
    !     ! type(ESMF_State)            ,intent(inout) :: aero_aci
    !     integer, intent(inout) :: aero_aci
    !     real, dimension (IM,JM,LM)  ,intent(in ) :: plo ! Pa
    !     real, dimension (IM,JM,0:LM),intent(in ) :: ple ! Pa
    !     real, dimension (IM,JM,LM)  ,intent(in ) :: q,t,omega,zlo, qlcn, qicn, qlls, qils
    !     real, dimension (IM,JM,0:LM),intent(in ) :: zle
    !     real, dimension (IM,JM)     ,intent(in ) :: FRLAND
    !     real, dimension (IM,JM)     ,intent(in ) :: sh, evap, kpbl     
    !     logical                     ,intent(in ) :: USE_AERO_BUFFER

    !     real, dimension (IM,JM,LM),intent(OUT) :: NACTL,NACTI, NWFA

    !     real(AER_PR), allocatable, dimension (:) :: sig0,rg,ni,bibar,nact 
    !     real(AER_PR), dimension (IM,JM)  :: naer_cb, zws
    !     real(AER_PR)                     :: wupdraft,tk,press,air_den,QC,QL,WC,BB,RAUX

    !     integer, dimension (IM,JM) :: kpbli

    !     real, dimension(:,:,:,:,:), allocatable :: buffer

    !     character(len=ESMF_MAXSTR)      :: aci_field_name
    !     real, pointer, dimension(:,:)   :: aci_ptr_2d
    !     real, pointer, dimension(:,:,:) :: aci_ptr_3d
    !     real, pointer, dimension(:,:,:) :: aci_num
    !     real, pointer, dimension(:,:,:) :: aci_dgn
    !     real, pointer, dimension(:,:,:) :: aci_sigma
    !     real, pointer, dimension(:,:,:) :: aci_density
    !     real, pointer, dimension(:,:,:) :: aci_hygroscopicity
    !     real, pointer, dimension(:,:,:) :: aci_f_dust
    !     real, pointer, dimension(:,:,:) :: aci_f_soot
    !     real, pointer, dimension(:,:,:) :: aci_f_organic
    !     character(len=ESMF_MAXSTR), allocatable, dimension(:) :: aero_aci_modes
    !     integer                         :: ACI_STATUS
    !     REAL                            :: aux1,aux2,aux3,hfs,hfl,nfaux

    !     integer :: n_modes
    !     REAL :: numbinit
    !     integer :: i,j,k,n,rc

    !     character(len=ESMF_MAXSTR)              :: IAm="Aer_Activation"
    !     integer                                 :: STATUS

    !     character(2) :: n_str
    !     character*100 :: dirName, rank_str
    !     integer :: rank, fileID

    !     do k = 1, LM
    !     do j = 1, JM
    !     do i = 1, IM
    !     AeroProps(i,j,k)%num = 0.0
    !     end do
    !     end do
    !     end do    
    !     NACTL    = 0.
    !     NACTI    = 0.

    !     kpbli = MAX(MIN(NINT(kpbl),LM-1),1)

    !     if (USE_AEROSOL_NN) then

    !     ! call ESMF_AttributeGet(aero_aci, name='number_of_aerosol_modes', value=n_modes, __RC__)

    !         open(newunit=fileID, file=trim(dirName) // "/n_modes_" // trim(rank_str) // ".in", &
    !             status='old', form="unformatted", action="read")
    !         read(fileID) n_modes
    !         close(fileID)

    !         if (n_modes > 0) then

    !     ! allocate( sig0(n_modes), __STAT__)
    !     ! allocate(   rg(n_modes), __STAT__)
    !     ! allocate(   ni(n_modes), __STAT__)
    !     ! allocate(bibar(n_modes), __STAT__)
    !     ! allocate( nact(n_modes), __STAT__)

    !     ! allocate(aero_aci_modes(n_modes), __STAT__)
    !     ! call ESMF_AttributeGet(aero_aci, name='aerosol_modes', itemcount=n_modes, valuelist=aero_aci_modes, __RC__)

    !     ! call ESMF_AttributeGet(aero_aci, name='air_pressure_for_aerosol_optics', value=aci_field_name, __RC__)
    !     ! if (aci_field_name /= '') then
    !     ! call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
    !     ! aci_ptr_3d = PLE
    !     ! end if

    !     ! call ESMF_AttributeGet(aero_aci, name='air_temperature', value=aci_field_name, __RC__)
    !     ! if (aci_field_name /= '') then
    !     ! call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
    !     ! aci_ptr_3d = T
    !     ! end if

    !     ! call ESMF_AttributeGet(aero_aci, name='fraction_of_land_type', value=aci_field_name, __RC__)
    !     ! if (aci_field_name /= '') then
    !     ! call MAPL_GetPointer(aero_aci, aci_ptr_2d, trim(aci_field_name), __RC__)
    !     ! aci_ptr_2d = FRLAND
    !     ! end if

    !             if (USE_AERO_BUFFER) then
    !             ! allocate(buffer(im,jm,lm,n_modes,8), __STAT__)
    !                 allocate(buffer(im,jm,lm,n_modes,8))
    !             end if

    !             allocate( sig0(n_modes))
    !             allocate(   rg(n_modes))
    !             allocate(   ni(n_modes))
    !             allocate(bibar(n_modes))
    !             allocate( nact(n_modes))

    !             allocate(aero_aci_modes(n_modes))

    !             ! *** Allocations done for standalone ***
    !             allocate(aci_num            (im, jm, lm))
    !             allocate(aci_dgn            (im, jm, lm))
    !             allocate(aci_sigma          (im, jm, lm))
    !             allocate(aci_density        (im, jm, lm))
    !             allocate(aci_hygroscopicity (im, jm, lm))
    !             allocate(aci_f_dust         (im, jm, lm))
    !             allocate(aci_f_soot         (im, jm, lm))
    !             allocate(aci_f_organic      (im, jm, lm))

    !             ACTIVATION_PROPERTIES: do n = 1, n_modes


    !                 if(n.ge.10) then
    !                     write(n_str,'(i2)') n
    !                 else
    !                     write(n_str,'(i1)') n
    !                 endif

    !     ! call ESMF_AttributeSet(aero_aci, name='aerosol_mode', value=trim(aero_aci_modes(n)), __RC__)

    !     ! ! execute the aerosol activation properties method 
    !     ! call ESMF_MethodExecute(aero_aci, label='aerosol_activation_properties', userRC=ACI_STATUS, RC=STATUS)
    !     ! VERIFY_(ACI_STATUS)
    !     ! VERIFY_(STATUS)

    !     ! ! copy out aerosol activation properties
    !     ! call ESMF_AttributeGet(aero_aci, name='aerosol_number_concentration', value=aci_field_name, __RC__)
    !     ! call MAPL_GetPointer(aero_aci, aci_num, trim(aci_field_name), __RC__)

    !     ! call ESMF_AttributeGet(aero_aci, name='aerosol_dry_size', value=aci_field_name, __RC__)
    !     ! call MAPL_GetPointer(aero_aci, aci_dgn, trim(aci_field_name), __RC__)

    !     ! call ESMF_AttributeGet(aero_aci, name='width_of_aerosol_mode', value=aci_field_name, __RC__)
    !     ! call MAPL_GetPointer(aero_aci, aci_sigma, trim(aci_field_name), __RC__)

    !     ! call ESMF_AttributeGet(aero_aci, name='aerosol_density', value=aci_field_name, __RC__)
    !     ! call MAPL_GetPointer(aero_aci, aci_density, trim(aci_field_name), __RC__)

    !     ! call ESMF_AttributeGet(aero_aci, name='aerosol_hygroscopicity', value=aci_field_name, __RC__)
    !     ! call MAPL_GetPointer(aero_aci, aci_hygroscopicity, trim(aci_field_name), __RC__)

    !     ! call ESMF_AttributeGet(aero_aci, name='fraction_of_dust_aerosol', value=aci_field_name, __RC__)
    !     ! call MAPL_GetPointer(aero_aci, aci_f_dust, trim(aci_field_name), __RC__)

    !     ! call ESMF_AttributeGet(aero_aci, name='fraction_of_soot_aerosol', value=aci_field_name, __RC__)
    !     ! call MAPL_GetPointer(aero_aci, aci_f_soot, trim(aci_field_name), __RC__)

    !     ! call ESMF_AttributeGet(aero_aci, name='fraction_of_organic_aerosol', value=aci_field_name, __RC__)
    !     ! call MAPL_GetPointer(aero_aci, aci_f_organic, trim(aci_field_name), __RC__)

    !                 open(newunit=fileID, file=trim(dirName) // "/aci_num_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
    !                     status='old', form="unformatted", action="read")
    !                 read(fileID) aci_num
    !                 close(fileID)

    !                 open(newunit=fileID, file=trim(dirName) // "/aci_dgn_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
    !                     status='old', form="unformatted", action="read")
    !                 read(fileID) aci_dgn
    !                 close(fileID)

    !                 open(newunit=fileID, file=trim(dirName) // "/aci_sigma_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
    !                     status='old', form="unformatted", action="read")
    !                 read(fileID) aci_sigma
    !                 close(fileID)

    !                 open(newunit=fileID, file=trim(dirName) // "/aci_density_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
    !                     status='old', form="unformatted", action="read")
    !                 read(fileID) aci_density
    !                 close(fileID)

    !                 open(newunit=fileID, file=trim(dirName) // "/aci_hygroscopicity_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
    !                     status='old', form="unformatted", action="read")
    !                 read(fileID) aci_hygroscopicity
    !                 close(fileID)

    !                 open(newunit=fileID, file=trim(dirName) // "/aci_f_dust_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
    !                     status='old', form="unformatted", action="read")
    !                 read(fileID) aci_f_dust
    !                 close(fileID)

    !                 open(newunit=fileID, file=trim(dirName) // "/aci_f_soot_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
    !                     status='old', form="unformatted", action="read")
    !                 read(fileID) aci_f_soot
    !                 close(fileID)

    !                 open(newunit=fileID, file=trim(dirName) // "/aci_f_organic_" // trim(n_str) // "_" // trim(rank_str) // ".in", &
    !                     status='old', form="unformatted", action="read")
    !                 read(fileID) aci_f_organic
    !                 close(fileID)

    !                 if (USE_AERO_BUFFER) then
    !                     buffer(:,:,:,n,1) = aci_num
    !                     buffer(:,:,:,n,2) = aci_dgn
    !                     buffer(:,:,:,n,3) = aci_sigma
    !                     buffer(:,:,:,n,4) = aci_hygroscopicity
    !                     buffer(:,:,:,n,5) = aci_density
    !                     buffer(:,:,:,n,6) = aci_f_dust
    !                     buffer(:,:,:,n,7) = aci_f_soot
    !                     buffer(:,:,:,n,8) = aci_f_organic
    !                 else
    !                     AeroProps(:,:,:)%num(n)   = aci_num
    !                     AeroProps(:,:,:)%dpg(n)   = aci_dgn
    !                     AeroProps(:,:,:)%sig(n)   = aci_sigma
    !                     AeroProps(:,:,:)%kap(n)   = aci_hygroscopicity
    !                     AeroProps(:,:,:)%den(n)   = aci_density
    !                     AeroProps(:,:,:)%fdust(n) = aci_f_dust
    !                     AeroProps(:,:,:)%fsoot(n) = aci_f_soot
    !                     AeroProps(:,:,:)%forg(n)  = aci_f_organic
    !                     AeroProps(:,:,:)%nmods    = n_modes                 ! no need of a 3D field: aero provider specific
    !                 end if

    !             end do ACTIVATION_PROPERTIES

    !             if (USE_AERO_BUFFER) then
    !                 do k = 1, LM
    !                     do j = 1, JM
    !                         do i = 1, IM
    !                             do n = 1, n_modes
    !                                 AeroProps(i,j,k)%num(n)   = buffer(i,j,k,n,1)
    !                                 AeroProps(i,j,k)%dpg(n)   = buffer(i,j,k,n,2)
    !                                 AeroProps(i,j,k)%sig(n)   = buffer(i,j,k,n,3)
    !                                 AeroProps(i,j,k)%kap(n)   = buffer(i,j,k,n,4)
    !                                 AeroProps(i,j,k)%den(n)   = buffer(i,j,k,n,5)
    !                                 AeroProps(i,j,k)%fdust(n) = buffer(i,j,k,n,6)
    !                                 AeroProps(i,j,k)%fsoot(n) = buffer(i,j,k,n,7)
    !                                 AeroProps(i,j,k)%forg(n)  = buffer(i,j,k,n,8)
    !                             end do
    !                             AeroProps(i,j,k)%nmods       = n_modes                 ! no need of a 3D field: aero provider specific
    !                         end do
    !                     end do
    !                 end do

    !                 deallocate(buffer)
    !             end if

    !             do k = 1, LM
    !                 do j = 1, JM
    !                     do i = 1, IM
    !                         nfaux =  0.0
    !                         do n = 1, n_modes
    !                             if (AeroProps(i,j,k)%kap(n) .gt. 0.4) then 
    !                                 nfaux =  nfaux + AeroProps(i,j,k)%num(n)
    !                             end if                           
    !                        end do !modes
    !                         NWFA(I, J, K)  =  nfaux
    !                     end do
    !                 end do 
    !             end do 

    !             deallocate(aero_aci_modes)

    !             !----- aerosol activation (single-moment uphysics)      
    !             do j = 1, JM
    !                 do i = 1, IM
    !                     aux1=  PLE(i,j,LM)/(287.04*(T(i,j,LM)*(1.+0.608*Q(i,j,LM)))) ! air_dens (kg m^-3)
    !                     hfs = -SH  (i,j) ! W m^-2
    !                     hfl = -EVAP(i,j) ! kg m^-2 s^-1
    !                     aux2= (hfs/MAPL_CP + 0.608*T(i,j,LM)*hfl)/aux1 ! buoyancy flux (h+le)
    !                     aux3=  ZLE(i,j,kpbli(i,j))           ! pbl height (m)
    !                     !-convective velocity scale W* (m/s)
    !                     ZWS(i,j) = max(0.,0.001-1.5*0.41*MAPL_GRAV*aux2*aux3/T(i,j,LM))
    !                     ZWS(i,j) = 1.2*ZWS(i,j)**0.3333 ! m/s           
    !                 enddo
    !             enddo

    !             !--- activated aerosol # concentration for liq/ice phases (units: m^-3)
    !             numbinit = 0.
    !             WC       = 0.
    !             BB       = 0.
    !             RAUX     = 0.

    !             !--- determing aerosol number concentration at cloud base
    !             DO j=1,JM
    !                 Do i=1,IM 
    !                     k            = kpbli(i,j)
    !                     naer_cb(i,j) = zero_par
    !                     tk           = T(i,j,k)              ! K
    !                     press        = plo(i,j,k)            ! Pa     
    !                     air_den      = press*28.8e-3/8.31/tk ! kg/m3
    !                     DO n=1,n_modes
    !                         if (AeroProps(i,j,k)%dpg(n) .ge. 0.5e-6) &
    !                         naer_cb(i,j)= naer_cb(i,j) + AeroProps(i,j,k)%num(n) * air_den         
    !                         naer_cb(i,j)= naer_cb(i,j)* 1.e+6  ! #/cm3
    !                         naer_cb(i,j)= max(real(1.0e-1,AER_PR),min(naer_cb(i,j),100.0))
    !                     ENDDO
    !                 ENDDO
    !             ENDDO!;ENDDO

    !             DO k=LM,1,-1
    !                 DO j=1,JM
    !                     Do i=1,IM

    !                         tk                 = T(i,j,k)                         ! K
    !                         press              = plo(i,j,k)                       ! Pa   
    !                         air_den            = press*28.8e-3/8.31/tk            ! kg/m3
    !                         qc                 = qicn(i,j,k)+qils(i,j,k)*1.e+3    ! g/kg
    !                         ql                 = qlcn(i,j,k)+qlls(i,j,k)*1.e+3    ! g/kg

    !                         IF( plo(i,j,k) > 34000.0) THEN 

    !                             wupdraft           = -9.81*air_den*omega(i,j,k)     ! m/s - grid-scale only              

    !                             !--in the boundary layer, add Wstar
    !                             if(k >= kpbli(i,j) .and. k < LM)  wupdraft = wupdraft+zws(i,j) 

    !                             IF(wupdraft > 0.1 .AND. wupdraft < 100.) THEN 

    !                                 ni   (1:n_modes)    =   max(AeroProps(i,j,k)%num(1:n_modes)*air_den,  zero_par)  ! unit: [m-3]
    !                                 rg   (1:n_modes)    =   max(AeroProps(i,j,k)%dpg(1:n_modes)*0.5*1.e6, zero_par)  ! unit: [um]
    !                                 sig0 (1:n_modes)    =       AeroProps(i,j,k)%sig(1:n_modes)                      ! unit: [um]
    !                                 bibar(1:n_modes)    =   max(AeroProps(i,j,k)%kap(1:n_modes),          zero_par)                 

    !                                 IF( tk >= 245.0) then   
    !                                     call GetActFrac(           n_modes    &
    !                                         ,      ni(1:n_modes)   &  
    !                                         ,      rg(1:n_modes)   & 
    !                                         ,    sig0(1:n_modes)   &  
    !                                         ,      tk              &
    !                                         ,   press              & 
    !                                         ,wupdraft              & 
    !                                         ,    nact(1:n_modes)   &
    !                                         ,   bibar(1:n_modes)   &
    !                                         )

    !                                     numbinit     = 0.
    !                                     NACTL(i,j,k) = 0.
    !                                     DO n=1,n_modes
    !                                         numbinit     = numbinit    + AeroProps(i,j,k)%num(n)*air_den
    !                                         NACTL(i,j,k) = NACTL(i,j,k)+ nact(n)
    !                                     ENDDO

    !                                     NACTL(i,j,k) = MIN(NACTL(i,j,k),0.99*numbinit)

    !                                 ENDIF ! tk>245
    !                             ENDIF   ! updraft > 0.1
    !                         ENDIF   ! plo > 34000.0


    !                         IF( tk <= 268.0) then
    !                             IF( (QC >= 0.5) .and. (QL >= 0.5)) then
    !                                 ! Number of activated IN following deMott (2010) [#/m3]  
    !                                 NACTI(i,j,k) = (1.e+3*ai*((273.16-tk)**bi) *  (naer_cb(i,j))**(ci*(273.16-tk)+di))  !#/m3
    !                             ELSE   !tk<243
    !                                 ! Number of activated IN following Wyser  

    !                                 WC    = air_den*QC  !kg/m3
    !                                 if (WC >= tiny(1.0)) then
    !                                     BB     =  -2. + log10(1000.*WC/50.)*(1.e-3*(273.15-tk)**1.5)
    !                                 else
    !                                     BB = -6.
    !                                 end if
    !                                 BB     = MIN((MAX(BB,-6.)),-2.)  

    !                                 RAUX   = 377.4 + 203.3 * BB+ 37.91 * BB **2 + 2.3696 * BB **3
    !                                 RAUX   = (betai + (gamai + deltai * RAUX**3)**0.5)**0.33333
    !                                 NACTI(i,j,k) = (3.* WC)/(4.*MAPL_PI*densic*(1.D-6*RAUX)**3)  !#/m3

    !                             ENDIF  !Mixed phase

    !                         ENDIF ! tk<=268
    !                         !
    !                         !
    !                         !-- fix limit for NACTL/NACTI
    !                         IF(NACTL(i,j,k) < NN_MIN) NACTL(i,j,k) = FRLAND(i,j)*NN_LAND + (1.0-FRLAND(i,j))*NN_OCEAN

    !                     ENDDO
    !                 ENDDO
    !             ENDDO

    !             deallocate(   rg)
    !             deallocate(   ni)
    !             deallocate(bibar)
    !             deallocate( nact)

    !         end if
    !     end if

    ! END SUBROUTINE Aer_Activation

!     subroutine GetActFrac(nmodes  & !nmodes                                     &
!         ,xnap                         & !ni              (1:nmodes)              & 
!         ,rg                           & !0.5d+00*dgn_dry (1:nmodes)              & 
!         ,sigmag                       & !sig0             (1:nmodes)              & 
!         ,tkelvin                      & !tk              (i,j,k)                     &
!         ,ptot                         & !pres             (i,j,k)                     & 
!         ,wupdraft                     & !wupdraft             (i,j,k)                     & 
!         ,nact                         & !nact             (i,j,k,1:nmodes)             &
!         ,bibar)
        
!         IMPLICIT NONE
  
!         ! arguments.
        
!         integer :: nmodes               !< number of modes [1]      
!         real(AER_PR) :: xnap(nmodes)         !< number concentration for each mode [#/m^3]
!         real(AER_PR) :: rg(nmodes)           !< geometric mean dry radius for each mode [um]
!         real(AER_PR) :: sigmag(nmodes)       !< geometric standard deviation for each mode [um]
!         real(AER_PR) :: tkelvin              !< absolute temperature [k]
!         real(AER_PR) :: ptot                 !< ambient pressure [pa]
!         real(AER_PR) :: wupdraft             !< updraft velocity [m/s]
!   !     real(AER_PR) :: ac(nmodes)           !< minimum dry radius for activation for each mode [um]
!   !     real(AER_PR) :: fracactn(nmodes)     !< activating fraction of number conc. for each mode [1]
!         real(AER_PR) :: nact(nmodes)         !< activating number concentration for each mode [#/m^3]
!         real(AER_PR) :: bibar(nmodes)        ! hygroscopicity parameter for each mode [1]
  
!         ! local variables. 
  
!         integer :: i, j                 ! loop counters       
  
!         !--------------------------------------------------------------------------------------------------------------
!         ! calculate the droplet activation parameters for each mode. 
!         !--------------------------------------------------------------------------------------------------------------
!         call ActFrac_Mat(nmodes,xnap,rg,sigmag,bibar,tkelvin,ptot,wupdraft,nact)
  
!     end subroutine GetActFrac

    ! subroutine ActFrac_Mat(nmodes,xnap,rg,sigmag,bibar,tkelvin,ptot,wupdraft,nact)

    !     IMPLICIT NONE
    
    !     ! Arguments.
        
    !     integer :: nmodes            !< number of modes [1]      
    !     real(AER_PR) :: xnap(nmodes)      !< number concentration for each mode [#/m^3]
    ! !     real(AER_PR) :: xmap(nmodes)      !< mass   concentration for each mode [ug/m^3]
    !     real(AER_PR) :: rg(nmodes)        !< geometric mean radius for each mode [um]
    !     real(AER_PR) :: sigmag(nmodes)    !< geometric standard deviation for each mode [um]
    !     real(AER_PR) :: bibar(nmodes)     !< hygroscopicity parameter for each mode [1]
    !     real(AER_PR) :: tkelvin           !< absolute temperature [k]
    !     real(AER_PR) :: ptot              !< ambient pressure [pa]
    !     real(AER_PR) :: wupdraft          !< updraft velocity [m/s]
    !     real(AER_PR) :: ac(nmodes)        !< minimum dry radius for activation for each mode [um]
    !     real(AER_PR) :: fracactn(nmodes)  !< activating fraction of number conc. for each mode [1]
    !     real(AER_PR) :: nact(nmodes)      !< activating number concentration for each mode [#/m^3]
    
    !     ! parameters.
        
    !     real(AER_PR), parameter :: pi            = 3.141592653589793d+00
    !     real(AER_PR), parameter :: twopi         = 2.0d+00 * pi
    !     real(AER_PR), parameter :: sqrt2         = 1.414213562d+00
    !     real(AER_PR), parameter :: threesqrt2by2 = 1.5d+00 * sqrt2
    
    !     real(AER_PR), parameter :: avgnum   = 6.0221367d+23       ! [1/mol]
    !     real(AER_PR), parameter :: rgasjmol = 8.31451d+00         ! [j/mol/k]
    !     real(AER_PR), parameter :: wmolmass = 18.01528d-03        ! molar mass of h2o     [kg/mol]
    !     real(AER_PR), parameter :: amolmass = 28.966d-03          ! molar mass of air     [kg/mol]
    !     real(AER_PR), parameter :: asmolmss = 132.1406d-03        ! molar mass of nh42so4 [kg/mol]
    !     real(AER_PR), parameter :: denh2o   = 1.00d+03            ! density of water [kg/m^3]
    !     real(AER_PR), parameter :: denamsul = 1.77d+03            ! density of pure ammonium sulfate [kg/m^3]
    !     real(AER_PR), parameter :: xnuamsul = 3.00d+00            ! # of ions formed when the salt is dissolved in water [1]
    !     real(AER_PR), parameter :: phiamsul = 1.000d+00           ! osmotic coefficient value in a-r 1998. [1] 
    !     real(AER_PR), parameter :: gravity  = 9.81d+00            ! grav. accel. at the earth's surface [m/s/s] 
    !     real(AER_PR), parameter :: heatvap  = 40.66d+03/wmolmass  ! latent heat of vap. for water and tnbp [j/kg] 
    !     real(AER_PR), parameter :: cpair    = 1006.0d+00          ! heat capacity of air [j/kg/k] 
    !     real(AER_PR), parameter :: t0dij    = 273.15d+00          ! reference temp. for dv [k] 
    !     real(AER_PR), parameter :: p0dij    = 101325.0d+00        ! reference pressure for dv [pa] 
    !     real(AER_PR), parameter :: dijh2o0  = 0.211d-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
    !     !----------------------------------------------------------------------------------------------------------------    
    !     ! real(AER_PR), parameter :: t0dij    = 283.15d+00          ! reference temp. for dv [k] 
    !     ! real(AER_PR), parameter :: p0dij    = 80000.0d+00         ! reference pressure for dv [pa] 
    !     ! real(AER_PR), parameter :: dijh2o0  = 0.300d-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
    !     !----------------------------------------------------------------------------------------------------------------
    !     real(AER_PR), parameter :: deltav   = 1.096d-07           ! vapor jump length [m]  
    !     real(AER_PR), parameter :: deltat   = 2.160d-07           ! thermal jump length [m]  
    !     real(AER_PR), parameter :: alphac   = 1.000d+00           ! condensation mass accommodation coefficient [1]  
    !     real(AER_PR), parameter :: alphat   = 0.960d+00           ! thermal accommodation coefficient [1]  
    
    !     ! local variables. 
    
    !     integer            :: i                              ! loop counter 
    !     real(AER_PR)            :: dv                             ! diffusion coefficient for water [m^2/s] 
    !     real(AER_PR)            :: dvprime                        ! modified diffusion coefficient for water [m^2/s] 
    !     real(AER_PR)            :: dumw, duma                     ! scratch variables [s/m] 
    !     real(AER_PR)            :: wpe                            ! saturation vapor pressure of water [pa]  
    !     real(AER_PR)            :: surten                         ! surface tension of air-water interface [j/m^2] 
    !     real(AER_PR)            :: xka                            ! thermal conductivity of air [j/m/s/k]  
    !     real(AER_PR)            :: xkaprime                       ! modified thermal conductivity of air [j/m/s/k]  
    !     real(AER_PR)            :: eta(nmodes)                    ! model parameter [1]  
    !     real(AER_PR)            :: zeta                           ! model parameter [1]  
    !     real(AER_PR)            :: xlogsigm(nmodes)               ! ln(sigmag) [1]   
    !     real(AER_PR)            :: a                              ! [m]
    !     real(AER_PR)            :: g                              ! [m^2/s]   
    !     real(AER_PR)            :: rdrp                           ! [m]   
    !     real(AER_PR)            :: f1                             ! [1]   
    !     real(AER_PR)            :: f2                             ! [1]
    !     real(AER_PR)            :: alpha                          ! [1/m]
    !     real(AER_PR)            :: gamma                          ! [m^3/kg]   
    !     real(AER_PR)            :: sm(nmodes)                     ! [1]   
    !     real(AER_PR)            :: dum                            ! [1/m]    
    !     real(AER_PR)            :: u                              ! argument to error function [1]
    !     real(AER_PR)            :: erf                            ! error function [1], but not declared in an f90 module 
    !     real(AER_PR)            :: smax                           ! maximum supersaturation [1]
    
    !     real(AER_PR)            :: aux1,aux2                           ! aux
    ! !----------------------------------------------------------------------------------------------------------------------
    ! !     rdrp is the radius value used in eqs.(17) & (18) and was adjusted to yield eta and zeta 
    ! !     values close to those given in a-z et al. 1998 figure 5. 
    ! !----------------------------------------------------------------------------------------------------------------------
    !     rdrp = 0.105d-06   ! [m] tuned to approximate the results in figures 1-5 in a-z et al. 1998.  
    ! !----------------------------------------------------------------------------------------------------------------------
    ! !     these variables are common to all modes and need only be computed once. 
    ! !----------------------------------------------------------------------------------------------------------------------
    !     dv = dijh2o0*(p0dij/ptot)*(tkelvin/t0dij)**1.94d+00                 ! [m^2/s] (p&k,2nd ed., p.503)
    !     surten = 76.10d-03 - 0.155d-03 * (tkelvin-273.15d+00)               ! [j/m^2] 
    !     wpe = exp( 77.34491296d+00 - 7235.424651d+00/tkelvin - 8.2d+00*log(tkelvin) + tkelvin*5.7113d-03 )  ! [pa] 
    !     dumw = sqrt(twopi*wmolmass/rgasjmol/tkelvin)                        ! [s/m] 
    !     dvprime = dv / ( (rdrp/(rdrp+deltav)) + (dv*dumw/(rdrp*alphac)) )   ! [m^2/s] - eq. (17) 
    !     xka = (5.69d+00+0.017d+00*(tkelvin-273.15d+00))*418.4d-05           ! [j/m/s/k] (0.0238 j/m/s/k at 273.15 k)
    !     duma = sqrt(twopi*amolmass/rgasjmol/tkelvin)                        ! [s/m]
    !     xkaprime = xka / ( ( rdrp/(rdrp+deltat) ) + ( xka*duma/(rdrp*alphat*denh2o*cpair) ) )   ! [j/m/s/k]
    !     g = 1.0d+00 / ( (denh2o*rgasjmol*tkelvin) / (wpe*dvprime*wmolmass) &
    !                     + ( (heatvap*denh2o) / (xkaprime*tkelvin) ) &
    !                     * ( (heatvap*wmolmass) / (rgasjmol*tkelvin) - 1.0d+00 ) )               ! [m^2/s]
    !     a = (2.0d+00*surten*wmolmass)/(denh2o*rgasjmol*tkelvin)                                 ! [m] 
    !     alpha = (gravity/(rgasjmol*tkelvin))*((wmolmass*heatvap)/(cpair*tkelvin) - amolmass)    ! [1/m] 
    !     gamma = (rgasjmol*tkelvin)/(wpe*wmolmass) &
    !             + (wmolmass*heatvap*heatvap)/(cpair*ptot*amolmass*tkelvin)                        ! [m^3/kg]
    !     dum = sqrt(alpha*wupdraft/g)                  ! [1/m] 
    !     zeta = 2.d+00*a*dum/3.d+00                    ! [1] 
    ! !----------------------------------------------------------------------------------------------------------------
    !     ! write(1,'(a27,4d15.5)')'surten,wpe,a            =',surten,wpe,a
    !     ! write(1,'(a27,4d15.5)')'xka,xkaprime,dv,dvprime =',xka,xkaprime,dv,dvprime
    !     ! write(1,'(a27,4d15.5)')'alpha,gamma,g, zeta     =',alpha,gamma,g,zeta
    ! !----------------------------------------------------------------------------------------------------------------------
    ! !     these variables must be computed for each mode. 
    ! !----------------------------------------------------------------------------------------------------------------------
    !     xlogsigm(:) = log(sigmag(:))                                                    ! [1] 
    !     smax = 0.0d+00                                                                  ! [1]
        
    !     do i=1, nmodes
        
    !         sm(i) = ( 2.0d+00/sqrt(bibar(i)) ) * ( a/(3.0*rg(i)) )**1.5d+00               ! [1] 
    !         eta(i) = dum**3 / (twopi*denh2o*gamma*xnap(i))                                ! [1] 
    
    !         !--------------------------------------------------------------------------------------------------------------
    !         ! write(1,'(a27,i4,4d15.5)')'i,eta(i),sm(i) =',i,eta(i),sm(i)
    !         !--------------------------------------------------------------------------------------------------------------
    !         f1 = 0.5d+00 * exp(2.50d+00 * xlogsigm(i)**2)                                 ! [1] 
    !         f2 = 1.0d+00 +     0.25d+00 * xlogsigm(i)                                     ! [1] 
    !         smax = smax + (   f1*(  zeta  / eta(i)              )**1.50d+00 &
    !                         + f2*(sm(i)**2/(eta(i)+3.0d+00*zeta))**0.75d+00 ) / sm(i)**2  ! [1] - eq. (6)
    !     enddo 
    !     smax = 1.0d+00 / sqrt(smax)                                                     ! [1]
        
        
    !     !aux1=0.0D0; aux2=0.0D0
    !     do i=1, nmodes
    
    !         ac(i)       = rg(i) * ( sm(i) / smax )**0.66666666666666667d+00               ! [um]
    
    !         u           = log(ac(i)/rg(i)) / ( sqrt2 * xlogsigm(i) )                      ! [1]
    !         fracactn(i) = 0.5d+00 * (1.0d+00 - erf(u))                                    ! [1]
    !         nact(i)     = fracactn(i) * xnap(i)                                           ! [#/m^3]
    !         !--------------------------------------------------------------------------------------------------------------
    !         !aux1=aux1+ xnap(i); aux2=aux2+ nact(i) 
            
    !         !if(fracactn(i) .gt. 0.9999999d+00 ) then
    !         !   write(*,*)i,ac(i),u,fracactn(i),xnap(i)
    !         !  print*,' xxx',i,ac(i),u,fracactn(i),xnap(i)
    !         ! stop
    !         ! endif
    !         !--------------------------------------------------------------------------------------------------------------
    !     enddo 
    
    !     return
    ! end subroutine ActFrac_Mat

    SUBROUTINE cup_cloud_limits(name,ierrc,ierr,cap_inc,cap_max,entr_rate              &
        ,heo_cup,heso_cup,qo_cup,qeso_cup,po,po_cup,z_cup,heo,hkbo &
        ,entr_rate_2d,hcot,k22,kbmax,klcl,kbcon,ktop            &
        ,use_excess,zqexec,ztexec,xland &
        ,itf,ktf,its,ite, kts,kte)

        IMPLICIT NONE
        ! character *(*), intent (in)         ::     name
        character *(10), intent (in)         ::     name

        integer ,intent (in   )             ::    &
            itf,ktf,its,ite, kts,kte,use_excess
        real, intent (in) :: entr_rate

        real,    dimension (its:ite,kts:kte)  ,intent (in   )     ::    &
            heo_cup,heso_cup,po,po_cup,z_cup,heo,qo_cup,qeso_cup
        real,    dimension (its:ite)          ,intent (in   )     ::    &
            cap_max,cap_inc,xland
        real,    dimension (its:ite)          ,intent (in   )     ::    &
            zqexec,ztexec
        integer, dimension (its:ite)          ,intent (in   )     ::    &
            kbmax
        integer, dimension (its:ite)          ,intent (inout)     ::    &
            kbcon,ierr,ktop,klcl,k22
        character*128                    ,intent (inout) :: ierrc(its:ite)
        real, dimension (its:ite)        ,intent (inout) :: hkbo
        real, dimension (its:ite,kts:kte),intent (inout) :: entr_rate_2d,hcot
        !  local variables in this routine
        integer                              ::   i,k,k1,k2,kfinalzu
        real                                 ::   depth_neg_buoy,plus,hetest,dz,x_add,dbythresh,frh
        real   , dimension (kts:kte)         ::   dby
        integer, dimension (its:ite)         ::   start_level
        integer,Parameter :: FIND_KTOP_OPTION = 1 !0=original, 1=new

!$acc data create(start_level, dby)
        !~ dbythresh=1.0!0.95  ! the range of this parameter is 0-1, higher => lower
        !~ ! overshoting (cheque AA0 calculation)
        !~ ! rainfall is too sensible this parameter
        !~ ! for now, keep =1.
!$acc kernels
        hcot = 0.0
        dby  = 0.0
        start_level = 0
!$acc end kernels
        !--- Check if is there any inversion layer between k22 and klcl
        !
!$acc parallel
!$acc loop gang
        DO i=its,itf
            IF(ierr(i) /= 0) cycle
!$acc loop vector
            do k=k22(i),klcl(i)
                hcot(i,k) = HKBO(I) ! assumed no entraiment between these layers
                !~ if(hcot(i,k) < heso_cup(i,k) )then
                !~ ierr(i)=102
                !~ ierrc(i)="there is an inversion layer between k22 and klcl"
                !~ print*,'inv=',k,k22(i),klcl(i);flush(6)
                !~ cycle
                !~ endif
            enddo
        ENDdo
!$acc end parallel
            !
            !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
            !
!$acc kernels
        loop0: DO i=its,itf
            !-default value
            kbcon(i)=KBMAX(i)+3

            IF(ierr(i) /= 0) cycle

            start_level(i)=klcl(i)

            loop1:  Do WHILE(ierr(i) == 0)

                kbcon(i)=start_level(i)
                do k=start_level(i)+1,KBMAX(i)+3
                    dz=z_cup(i,k)-z_cup(i,k-1)
                    hcot(i,k)= ( (1.-0.5*entr_rate_2d(i,k-1)*dz)*hcot(i,k-1)     &
                                + entr_rate_2d(i,k-1)*dz *heo (i,k-1)   ) / &
                            (1.+0.5*entr_rate_2d(i,k-1)*dz)
                enddo

                !~ IF(NAME == 'shallow') then
                !~ loopsh:              do k=start_level(i)+1,KBMAX(i)+1
                !~ if(hcot(i,k) >= HESO_cup(i,k)) then
                !~ kbcon(i) = k
                !~ exit loopsh
                !~ endif
                !~ enddo loopsh

                !~ ELSE

                loop2:        DO while (hcot(i,kbcon(i)) < HESO_cup(i,kbcon(i)))                    
                    kbcon(i)=kbcon(i)+1
                    if(kbcon(i).gt.kbmax(i)+2) then
                        ierr(i)=3
                        ierrc(i)="could not find reasonable kbcon in cup_kbcon"
                        Exit loop2
                    endif
                !print*,"kbcon=",kbcon(i);flush(6)
                enddo loop2
                !~ ENDIF

                IF(ierr(i) /= 0) cycle loop0

                !---     cloud base pressure and max moist static energy pressure
                !---     i.e., the depth (in mb) of the layer of negative buoyancy
                depth_neg_buoy = - (po_cup(i,kbcon(i))-po_cup(i,start_level(i)))

                !- test if the air parcel has enough energy to reach the positive buoyant region
                if(cap_max(i) > depth_neg_buoy) cycle loop0

                !- if am here -> kbcon not found for air parcels from k22 level
                k22(i)=k22(i)+1
                !- get new hkbo
                x_add = float(use_excess)*(xlv*zqexec(i)+cp*ztexec(i))
                call get_cloud_bc(name,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup(i,kts:kte),hkbo(i),k22(i),x_add)
                !
                start_level(i)=start_level(i)+1
                !
                hcot(i,start_level(i))=hkbo (i)
            ENDDO loop1
        ENDDO loop0
!$acc end kernels
        !
        if(NAME /= 'shallow') return
        !
        !--- DETERMINE THE LEVEL OF NEUTRAL BUOYANCY - KTOP
        !
!$acc kernels
        DO i=its,itf
            kfinalzu=ktf-1
            ktop(i) =kfinalzu
            IF(ierr(i) /= 0) cycle
            !~ dby(:)=0.0

            start_level(i)=kbcon(i)

            !print*,'hco1=',start_level(I),kbcon(i),hcot(i,start_level(i))/heso_cup(i,start_level(i))

            do k=start_level(i)+1,ktf-1
                dz=z_cup(i,k)-z_cup(i,k-1)

                hcot(i,k)=( (1.-0.5*entr_rate_2d(i,k-1)*dz)*hcot(i,k-1)    &
                        +entr_rate_2d(i,k-1)*dz *heo (i,k-1) )/ &
                        (1.+0.5*entr_rate_2d(i,k-1)*dz)
                !~ dby(k)=dby(k-1)+(hcot(i,k)-heso_cup(i,k))*dz
                !print*,'hco2=',k,hcot(i,k)/heso_cup(i,k),dby(k),entr_rate_2d(i,k-1)

            enddo
            !~ if(FIND_KTOP_OPTION==0) then
            !~ do k=maxloc(dby(:),1),ktf-1
            !~ print*,'hco30=',k,dby(k),dbythresh*maxval(dby)

            !~ if(dby(k).lt.dbythresh*maxval(dby))then
            !~ kfinalzu = k - 1
            !~ ktop(i)  = kfinalzu
            !~ !print*,'hco4=',k,kfinalzu,ktop(i),kbcon(i)+1;flush(6)
            !~ go to 412
            !~ endif
            !~ enddo
            !~ 412    continue
            !~ else
            do k=start_level(i)+1,ktf-1
                !print*,'hco31=',k,hkbo(i),hcot(i,k),heso_cup(i,k)

                if(hcot(i,k) < heso_cup(i,k) )then
                    kfinalzu = k - 1
                    ktop(i)  = kfinalzu
                    !print*,'hco40=',k,kfinalzu,ktop(i),kbcon(i);flush(6)
                    exit
                endif
            enddo
            !~ endif
            !~ print*,'hco5=',k,kfinalzu,ktop(i),kbcon(i);flush(6)

            if(kfinalzu.le.kbcon(i)+1) ierr(i)=41
            !print*,"3",ierr(i),klcl(i),kbcon(i),ktop(i);flush(6)

        ENDDO
!$acc end kernels
!$acc end data
    END SUBROUTINE cup_cloud_limits

    SUBROUTINE get_inversion_layers(cumulus,ierr,psur,po_cup,to_cup,zo_cup,k_inv_layers,&
        dtempdz,itf,ktf,its,ite, kts,kte)

        IMPLICIT NONE
        INTEGER,                              INTENT (IN )  :: itf,ktf,its,ite,kts,kte
        ! CHARACTER *(*),                       INTENT (IN )  :: cumulus
        CHARACTER *(10),                       INTENT (IN )  :: cumulus
        INTEGER, DIMENSION (its:ite),         INTENT (INOUT):: ierr
        REAL,    DIMENSION (its:ite),         INTENT (IN )  :: psur
        REAL,    DIMENSION (its:ite,kts:kte), INTENT (IN )  :: po_cup,to_cup,zo_cup
        REAL,    DIMENSION (its:ite,kts:kte), INTENT (OUT)  :: dtempdz
        INTEGER, DIMENSION (its:ite,kts:kte), INTENT (OUT)  :: k_inv_layers
        REAL :: dzm,delp, first_deriv(kts:kte),sec_deriv(kts:kte),distance(kts:kte)
        INTEGER:: i,k,ilev,kk,k1,ix,k800,k550,ist
        INTEGER, PARAMETER :: extralayer = 0 !- makes plume top higher
        INTEGER, DIMENSION (its:ite,kts:kte)  :: local_k_inv_layers
        !
        !-initialize k_inv_layers as 1 (non-existent layer)_
        k_inv_layers= 1 !integer 
        dtempdz     = 0.0
        first_deriv = 0.0
        sec_deriv   = 0.0
        distance    = 0.0
        local_k_inv_layers=1
        ist=3
!$acc kernels
        DO i = its,itf
            if(ierr(i) /= 0) cycle
            !- displacement from local surface pressure level
            delp=1000.-psur(i)

            !- 2nd method
            ! DO k = kts+1,ktf-2
            !dtempdz(i,k)=   ( deriv3(zo_cup(i,k), zo_cup(i,kts:ktf), to_cup(i,kts:ktf), ktf-kts+1, 1,ierr(i)))
            !!! sec_deriv(k)=abs( deriv3(zo_cup(i,k), zo_cup(i,kts:ktf), to_cup(i,kts:ktf), ktf-kts+1, 2))
            !print*,"2=",k,z_cup(i,k),dtempdz(i,k),
            ! ENDDO
            ! if(ierr(i) /= 0) cycle

            !-1st method
            !-  get the 1st derivative
            DO k = kts+ist,ktf-ist
                first_deriv(k)  = (to_cup(i,k+1)-to_cup(i,k-1))/(zo_cup(i,k+1)-zo_cup(i,k-1))
            ENDDO
            first_deriv(kts      :kts+ist-1)  =first_deriv(kts+ist)
            first_deriv(ktf-ist+1:kte      )  =first_deriv(ktf-ist)

            dtempdz  (i,:)  = first_deriv(:)

            !-  get the abs of the 2nd derivative
            DO k = kts+ist+1,ktf-ist-1
                sec_deriv(k)= abs((first_deriv(k+1)-first_deriv(k-1))/(zo_cup(i,k+1)-zo_cup(i,k-1)))
            ENDDO
            sec_deriv(kts    :kts+ist)=sec_deriv(kts+ist+1)
            sec_deriv(ktf-ist:kte    )=sec_deriv(ktf-ist-1)

            ix=1
            DO kk=kts+ist+2,ktf-ist-2
                if(sec_deriv(kk) < sec_deriv(kk+1) .and. sec_deriv(kk) < sec_deriv(kk-1)) then
                    local_k_inv_layers(i,ix)=kk 
                    ix  =ix+1
                endif
            ENDDO

            !- 2nd criteria
            DO k=kts+ist+2,ktf-ist-2
                kk=local_k_inv_layers(i,k)
                if(kk == 1) cycle
                if( dtempdz(i,kk) < dtempdz(i,kk-1) .and. dtempdz(i,kk) < dtempdz(i,kk+1) ) then ! the layer is not a local maximum
                    local_k_inv_layers(i,k) = 1
                endif
            ENDDO

        ENDDO


        !- find the locations of inversions around 800 and 550 hPa
        DO i = its,itf
            !----------------
            !k_inv_layers(i,mid)=1
            !----------------
            if(ierr(i) /= 0) cycle
            !- displacement from local surface pressure level
            delp=1000.-psur(i)
            !----------------
            !k_inv_layers(i,mid)=21
            !cycle
            !----------------
            ! IF( trim(cumulus)=='shallow') THEN 
            IF( cumulus(1:7)=='shallow') THEN
                !- now find the closest layers of 800 and 550 hPa.
                !- this is for shallow convection k800
                DO k=kts,ktf
                    distance(k)=abs(po_cup(i,local_k_inv_layers(i,k))-(800.-delp))
                ENDDO
                k800=minloc(abs(distance(kts:ktf)),1)

                if( k800 <= kts .or. k800 >= ktf - 4) then
                    k_inv_layers(i,shal)= 1
                    ierr(i)=8
                else
                    !-save k800 in the k_inv_layers array
                    k_inv_layers(i,shal)=local_k_inv_layers(i,k800) +extralayer
                endif
                if(  k_inv_layers(i,shal) <= kts .or. k_inv_layers(i,shal) >= ktf-4) then 
                    !print*,"SHAL_k_inv_layers=",k_inv_layers(i,shal),ierr(i)
                    ierr(i)=11
                endif

            ! ELSEIF( trim(cumulus)=='mid') THEN 
            ELSEIF( cumulus(1:3)=='mid') THEN 
                !- this is for mid/congestus convection k500
                DO k=kts,ktf
                    distance(k)=abs(po_cup(i,local_k_inv_layers(i,k))-(550.-delp))
                ENDDO
                k550=minloc(abs(distance(kts:ktf)),1)

                if( k550 <= kts .or. k550 >= ktf - 4) then
                    k_inv_layers(i,mid) = 1
                    ierr(i)=8
                else
                    !-save k550 in the k_inv_layers array
                    k_inv_layers(i,mid )=local_k_inv_layers(i,k550) +extralayer
                endif
                if(  k_inv_layers(i,mid) <= kts .or. k_inv_layers(i,mid) >= ktf-4) then 
                    !print*,"MID_k_inv_layers=",k_inv_layers(i,MID),ierr(i)
                    ierr(i)=12
                endif
            ELSE
                k_inv_layers(i,:)=1
                ierr(i)=88
            ENDIF

            !-set the rest to undef (not in use)
            !where(k_inv_layers(i,:) =/ 1) k_inv_layers(i,:)=1

            !k_inv_layers(i,4:kte)= 1
            !k_inv_layers(i,kts  )= 1

            !~ write(13,100) k_inv_layers(i,1:maxloc(k_inv_layers(i,:),1))
            !~ write(13,101) po_cup(i,k_inv_layers(i,1:7))
            !~ write(14,102) 'shall',k_inv_layers(i,shal),po_cup(i,k_inv_layers(i,shal))
            !~ write(13,102) 'mid  ',k_inv_layers(i,mid),po_cup(i,k_inv_layers(i,mid))
            !~ 100 format(1x,8i8)
            !~ 101 format(1x,8F10.3)
            !~ 102 format(1x,A7,i8,F10.3)

        ENDDO
!$acc end kernels
        !-only for debugging
        !~ do i=its,itf
        !~ if(ierr(i) /= 0) cycle
        !~ if((k_inv_layers(i,mid)) .ge. kte .or. (k_inv_layers(i,mid)) .le. kts) then
        !~ print*,"wrong inv layr",kts,kte,k_inv_layers(i,mid),mid
        !~ !stop "INV LAYER"
        !~ endif
        !~ ENDDO

! contains
! real function deriv3(xx, xi, yi, ni, m,ierr)
! !====================================================================
! ! Evaluate first- or second-order derivatives
! ! using three-point Lagrange interpolation
! ! written by: Alex Godunov (October 2009)
! !--------------------------------------------------------------------
! ! input ...
! ! xx    - the abscissa at which the interpolation is to be evaluated
! ! xi()  - the arrays of data abscissas
! ! yi()  - the arrays of data ordinates
! ! ni - size of the arrays xi() and yi()
! ! m  - order of a derivative (1 or 2)
! ! output ...
! ! deriv3  - interpolated value
! !============================================================================*/

! implicit none
! integer, parameter :: n=3
! real   , intent(in):: xx
! integer, intent(in):: ni, m
! real   , intent(in) :: xi(ni), yi(ni)
! real:: x(n), f(n)
! integer i, j, k, ix
! integer, intent(inout) :: ierr

! ! exit if too high-order derivative was needed,
! if (m > 2) then
! deriv3 = 0.0
! return
! end if

! ! if x is ouside the xi(1)-xi(ni) interval set deriv3=0.0
! if (xx < xi(1) .or. xx > xi(ni)) then
! deriv3 = 0.0
! ierr=8
! !stop "problem with 2nd derivative-deriv3 routine"
! return
! endif

! ! a binary (bisectional) search to find i so that xi(i-1) < x < xi(i)
! i = 1
! j = ni
! do while (j > i+1)
! k = (i+j)/2
! if (xx < xi(k)) then
! j = k
! else
! i = k
! endif
! enddo

! ! shift i that will correspond to n-th order of interpolation
! ! the search point will be in the middle in x_i, x_i+1, x_i+2 ...
! i = i + 1 - n/2

! ! check boundaries: if i is ouside of the range [1, ... n] -> shift i
! if (i < 1) i=1
! if (i + n > ni) i=ni-n+1

! !  old output to test i
! !  write(*,100) xx, i
! !  100 format (f10.5, I5)

! ! just wanted to use index i
! ix = i

! ! initialization of f(n) and x(n)
! do i=1,n
! f(i) = yi(ix+i-1)
! x(i) = xi(ix+i-1)
! end do

! ! calculate the first-order derivative using Lagrange interpolation
! if (m == 1) then
! deriv3 =          (2.0*xx - (x(2)+x(3)))*f(1)/((x(1)-x(2))*(x(1)-x(3)))
! deriv3 = deriv3 + (2.0*xx - (x(1)+x(3)))*f(2)/((x(2)-x(1))*(x(2)-x(3)))
! deriv3 = deriv3 + (2.0*xx - (x(1)+x(2)))*f(3)/((x(3)-x(1))*(x(3)-x(2)))
! ! calculate the second-order derivative using Lagrange interpolation
! else
! deriv3 =          2.0*f(1)/((x(1)-x(2))*(x(1)-x(3)))
! deriv3 = deriv3 + 2.0*f(2)/((x(2)-x(1))*(x(2)-x(3)))
! deriv3 = deriv3 + 2.0*f(3)/((x(3)-x(1))*(x(3)-x(2)))
! endif
! end function deriv3
    END SUBROUTINE get_inversion_layers

    SUBROUTINE get_zu_zd_pdf( &
        cumulus, &
        draft,&
        ierr,&
        kb,&
        kt,&
        zu,&
        kts,&
        kte,&
        ktf,&
        kpbli,&
        k22,&
        kbcon,&
        klcl,&
        po_cup,&
        psur,&
        xland)
        !$acc routine seq
        implicit none
        integer, intent(in   ) :: kts,kte,ktf,kpbli,k22,kbcon,kt,kb,klcl
        integer, intent(inout) :: ierr
        real   , intent(in   ) :: po_cup(kts:kte)
        real   , intent(in   ) :: psur,xland
        real   , intent(inout) :: zu(kts:kte)
        character*(*), intent(in) ::draft,cumulus
        !- local var
        integer :: kk,add,i,nrec=0,k,kb_adj,kpbli_adj,level_max_zu,ktarget
        real :: zumax,ztop_adj,a2,beta, alpha,kratio,tunning,FZU,krmax,dzudk
        real :: zuh(kts:kte),zul(kts:kte),  pmaxzu ! pressure height of max zu for deep
        real,   parameter :: px =45./120. ! px sets the pressure level of max zu. its range is from 1 to 120.
        real,   parameter :: px2=45./120. ! px sets the pressure level of max zu. its range is from 1 to 120.
        integer:: itest                ! 5=gamma+beta, 4=gamma, 1=beta
      
        integer:: minzu,maxzul,maxzuh
        logical, parameter :: do_smooth = .TRUE. 
      
        !-------- gama pdf
        real, parameter :: beta_deep=1.25,g_beta_deep=0.8974707
        INTEGER :: k1
        real :: lev_start,g_alpha2,g_a,y1,x1,g_b,a,b,alpha2,csum,zubeg,wgty
        real, dimension(30) :: x_alpha,g_alpha
        data  (x_alpha(k),k=1,30)/                                    & 
                      3.699999,3.699999,3.699999,3.699999,            &
                        3.024999,2.559999,2.249999,2.028571,1.862500, &
                        1.733333,1.630000,1.545454,1.475000,1.415385, &
                        1.364286,1.320000,1.281250,1.247059,1.216667, &
                        1.189474,1.165000,1.142857,1.122727,1.104348, &
                        1.087500,1.075000,1.075000,1.075000,1.075000, &
                      1.075000/
        data (g_alpha(k),k=1,30)/                                         &
                      4.1706450,4.1706450,4.1706450,4.1706450,            &
                        2.0469250,1.3878370,1.1330030,1.012418,0.9494680, &
                        0.9153771,0.8972442,0.8885444,0.8856795,0.8865333,&
                        0.8897996,0.8946404,0.9005030,0.9070138,0.9139161,&
                        0.9210315,0.9282347,0.9354376,0.9425780,0.9496124,&
                        0.9565111,0.9619183,0.9619183,0.9619183,0.9619183,&
                      0.9619183/
        ! -------- gama pdf
      
        !-- fill zu with zeros
        itest=-999
        zu =0.0
        zuh=0.0
        zul=0.0
        if(draft(1:7) == "deep_up" .and. xland >  0.90) itest=11 !ocean
        if(draft(1:7) == "deep_up" .and. xland <= 0.90) itest=12 !land
        if(draft(1:6) == "mid_up" ) itest= 5
      
        ! ---------------------------------------------------------
        IF(itest==5 .and. draft(1:6) == "mid_up" ) then
      
            !--- part 1 GAMMA format
            csum =0.
            zubeg=0.
            lev_start=min(.9,.1+csum*.013)
            kb_adj=max(kb,2)
            kb_adj=min(kb_adj,kt-1)
            if(kb_adj==kt) stop "kb_adj==kt"
      
            tunning = 0.30
            alpha2  = (tunning*(beta_deep -2.)+1.)/(1.-tunning)
      
            do k=27,3,-1
                  if(x_alpha(k) >= alpha2)exit
            enddo
            k1=k+1
            if(x_alpha(k1) .ne. x_alpha(k1-1))then
              a=x_alpha(k1)-x_alpha(k1-1)
              b=x_alpha(k1-1)*(k1) -(k1-1)*x_alpha(k1)
              x1= (alpha2-b)/a
              y1=a*x1+b
              g_a=g_alpha(k1)-g_alpha(k1-1)
              g_b=g_alpha(k1-1)*k1 - (k1-1)*g_alpha(k1)
              g_alpha2=g_a*x1+g_b
            else
              g_alpha2=g_alpha(k1)
            endif
      
            fzu = gamma(alpha2 + beta_deep)/(g_alpha2*g_beta_deep)
            fzu=0.01*fzu
            do k=kb_adj,min(kte,kt)
               kratio= (po_cup(k)-po_cup(kb_adj))/(po_cup(kt)-po_cup(kb_adj))
               zu(k) = zubeg+FZU*kratio**(alpha2-1.0) * (1.0-kratio)**(beta_deep-1.0)
              
            enddo
            !- normalize ZU
            zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ (1.e-12+maxval(zu(kts:min(kte,kt+1))))
      
            !--- part 2: BETA format
            pmaxzu=psur-px*(psur-po_cup(kt))
            ! kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
            kb_adj=minloc_acc(abs(po_cup(kts:kt)-pmaxzu))
            kb_adj=max(kb,kb_adj)
            kb_adj=min(kb_adj,kt)
            !beta=4.  !=> must be larger than 1
                      !=> higher makes the profile sharper
                      !=> around the maximum zu
            !- 2nd approach for beta and alpha parameters
            !- the tunning parameter must be between 0.5 (low  level max zu)
            !-                                   and 1.5 (high level max zu)
            !tunning= 1.0
            tunning = 0.6 !-- clo-X and tune0.6 experiment.
            !
            beta    = 2.0/tunning
            alpha   = tunning*beta
            !
            !- this alpha constrains the location of the maximun ZU to be at
            !- "kb_adj" vertical level
            alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
            !
            ! imposing zu(ktop) = 0
            do k=klcl-1,min(kte,kt)
                kratio= float(k)/float(kt+1)
                zuh(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
            enddo
            !- normalize ZU
            zuh(kts:min(kte,kt+1))= zuh(kts:min(kte,kt+1))/ (1.e-12+maxval(zuh(kts:min(kte,kt+1))))
      
            !--- part 3: BETA format from sfc to max zuh, then GAMMA format
            ! Do k=kts,max(kts,maxloc(zuh(:),1)-2)
            Do k=kts,max(kts,maxloc_acc(zuh(:))-2)
                zu(k)=zuh(k)
            ENDDO
            ! Do k=max(kts,maxloc(zuh(:),1)-1),min(maxloc(zuh(:),1)+1,kt)
            Do k=max(kts,maxloc_acc(zuh(:))-1),min(maxloc_acc(zuh(:))+1,kt)
                zu(k)=0.5*(zu(k)+zuh(k))
            ENDDO
      
            !-- special treatment below k22/klcl
            DO k=klcl,kts+1,-1
              zu(k)=zu(k+1)*0.5
            enddo
            !-- smooth section
            IF(do_smooth) then
                !--from surface
                zul(kts+1)=zu(kts+1)*0.25
                ! do k=kts+2,maxloc(zu,1)
                do k=kts+2,maxloc_acc(zu)
                   zul(k)=(zu(k-1)+zu(k))*0.5
                enddo
                ! do k=kts+1,maxloc(zu,1)
                do k=kts+1,maxloc_acc(zu)
                   zu(k)=(zul(k)+zu(k))*0.5
                enddo
      
                !--from ktop
                zul(kt-1)=zu(kt-1)*0.1
                !print*,"ZUMD=",kt,zu(kt),zul(kt)
                ! do k=kt-2,kt-min(maxloc(zu,1),5),-1
                do k=kt-2,kt-min(maxloc_acc(zu),5),-1
                   zul(k)=(zul(k+1)+zu(k))*0.5
                enddo
                wgty=0.
                ! do k=kt,kt-min(maxloc(zu,1),5),-1
                do k=kt,kt-min(maxloc_acc(zu),5),-1
                    ! wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
                    wgty=wgty+1./(float(min(maxloc_acc(zu),5))+1)
                    zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
                   !print*,"zuMD=",k,zu(k),zul(k),(zul(k)+zu(k))*0.5,min(maxloc(zu,1),5),wgty
                enddo
            ENDIF
            zu(kts)=0.
        !---------------------------------------------------------
        ELSEIF(itest==12 .and. draft(1:7) == "deep_up") then
            !- kb cannot be at 1st level
            !kb_adj=max(kb,2)
            pmaxzu=psur-px*(psur-po_cup(kt))
            ! kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
            kb_adj=minloc_acc(abs(po_cup(kts:kt)-pmaxzu))
            kb_adj=max(kb,kb_adj)
            kb_adj=min(kb_adj,kt)
      
            !beta=4.  !=> must be larger than 1
                      !=> higher makes the profile sharper
                      !=> around the maximum zu
            !- 2nd approach for beta and alpha parameters
            !- the tunning parameter must be between 0.5 (low  level max zu)
            !-                                   and 1.5 (high level max zu)
             tunning = 1.15  !NEW
            !tunning = 1.   !B7p1
            !tunning = 0.6  !P9P6
            beta    = 2.0/tunning
            alpha   = tunning*beta
            !
            !- this alpha constrains the location of the maximun ZU to be at
            !- "kb_adj" vertical level
            alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
            !
            ! imposing zu(ktop) = 0
            do k=klcl-1,min(kte,kt)
                kratio= float(k)/float(kt+1)
                zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
            enddo
      
            !-- zu = zero below k22
            !zu(1:max(k22-1,1)) = 0.0
            !-- special treatment below k22/klcl
            DO k=klcl,kts+1,-1
                zu(k)=zu(k+1)*0.5
            enddo
           !-- smooth section
           IF(do_smooth) then
                !--from surface
                zul(kts+1)=zu(kts+1)*0.25
                ! do k=kts+2,maxloc(zu,1)
                do k=kts+2,maxloc_acc(zu)
                    zul(k)=(zu(k-1)+zu(k))*0.5
                enddo
                ! do k=kts+1,maxloc(zu,1)
                do k=kts+1,maxloc_acc(zu)
                    zu(k)=(zul(k)+zu(k))*0.5
                enddo
      
                !--from ktop
                zul(kt)=zu(kt)*0.1
                ! do k=kt-1,kt-min(maxloc(zu,1),5),-1
                do k=kt-1,kt-min(maxloc_acc(zu),5),-1
                    zul(k)=(zul(k+1)+zu(k))*0.5
                enddo
                wgty=0.0
                ! do k=kt,kt-min(maxloc(zu,1),5),-1
                do k=kt,kt-min(maxloc_acc(zu),5),-1
                    ! wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
                    wgty=wgty+1./(float(min(maxloc_acc(zu),5))+1)
                    zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
                enddo
            ENDIF
      
            zu(kts)=0.
        !---------------------------------------------------------
        ELSEIF(itest==11 .and. draft(1:7) == "deep_up") then
            !pmaxzu=psur-px*(psur-po_cup(kt))
            pmaxzu=850.
            ! kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)!;print*,"1=",kb_adj,po_cup(kb_adj)
            kb_adj=minloc_acc(abs(po_cup(kts:kt)-pmaxzu))
            kb_adj=max(kb,kb_adj)
            kb_adj=min(kb_adj,kt)
            !tunning = 1.  !-- P10
            tunning = 0.6 !-- clo-X and tune0.6 experiment., B7p1
      
            beta    = 2.0/tunning
            !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
            alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
            do k=klcl-1,min(kte,kt)
                kratio= float(k)/float(kt+1)
                zul(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
                !print*,"1",k,zul(k),kb_adj,pmaxzu
            enddo
            zul(kts:min(kte,kt))= zul(kts:min(kte,kt))/ (1.e-9+maxval(zul(kts:min(kte,kt)),1))
      
            !-----------
            pmaxzu=po_cup(kt)+150.
            !pmaxzu=200.
            ! kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)!;print*,"2=",kb_adj,po_cup(kb_adj),tune
            kb_adj=minloc_acc(abs(po_cup(kts:kt)-pmaxzu))
            kb_adj=max(kb,kb_adj)
            kb_adj=min(kb_adj,kt)
            !tunning = 1.  !P10, B7p1
            !tunning = 0.6 !P9P6
            tunning = 0.8 !-- NEW
      
            beta    = 1.0/tunning
      
            !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
            alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
            do k=klcl-1,min(kte,kt)
                kratio= float(k)/float(kt+1)
                zuh(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
            enddo
            zuh(kts:min(kte,kt))= zuh(kts:min(kte,kt))/ (1.e-9+maxval(zuh(kts:min(kte,kt)),1))
      
            !increasing contribuition of zuh => more heating at upper levels
            !zu(:)=0.7*zul(:)+ 0.3*zuh(:) !P10, B7p1
            !zu(:)=0.9*zul(:)+ 0.1*zuh(:) !P9P6
            zu(:)=0.65*zul(:)+ 0.35*zuh(:) !NEW
      
            !-- special treatment below k22/klcl
            DO k=klcl,kts+1,-1
                zu(k)=zu(k+1)*0.5
            enddo
            !-- smooth section
            IF(do_smooth) then
                !--from surface
                zul(kts+1)=zu(kts+1)*0.25
                ! do k=kts+2,maxloc(zu,1)
                do k=kts+2,maxloc_acc(zu)
                   zul(k)=(zu(k-1)+zu(k))*0.5
                enddo
                ! do k=kts+1,maxloc(zu,1)
                do k=kts+1,maxloc_acc(zu)
                   zu(k)=(zul(k)+zu(k))*0.5
                enddo
      
                !--from ktop
                zul(kt)=zu(kt)*0.1
                ! do k=kt-1,kt-min(maxloc(zu,1),5),-1
                do k=kt-1,kt-min(maxloc_acc(zu),5),-1
                   zul(k)=(zul(k+1)+zu(k))*0.5
                enddo  
            
                wgty=0.
                ! do k=kt,kt-min(maxloc(zu,1),5),-1
                do k=kt,kt-min(maxloc_acc(zu),5),-1
                    ! wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
                    wgty=wgty+1./(float(min(maxloc_acc(zu),5))+1)
                    zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
                    !print*,"zu=",k,zu(k),zul(k),(zul(k)+zu(k))*0.5,min(maxloc(zu,1),5),wgty
                enddo
            ENDIF
            zu(kts)=0.
        !---------------------------------------------------------
        ELSEIF(draft(1:10) == "shallow_up") then
            !kb_adj=kb ! level where mass flux starts
            kb_adj=klcl-1!kbcon-1! level where mass flux starts
            kpbli_adj=kpbli
            if(kpbli_adj < kb_adj .or. kpbli_adj >= kt ) then
                kpbli_adj = kb_adj + 1
            endif
      
            !- location of the maximum Zu: 2 levels above
            level_max_zu=1
            krmax=float(kpbli_adj-kb_adj+1+level_max_zu)/float(kt+1)
            krmax=min(krmax,0.99)
            !
            beta= 4.
            !- this alpha imposes the maximum zu at kpbli
            alpha=1.+krmax*(beta-1.)/(1.-krmax)
            alpha=min(6.,alpha)
      
            !- to check if dZu/dk = 0 at k=kpbli_adj
            !kratio=krmax
            !dzudk=(alpha-1.)*(kratio**(alpha-2.)) * (1.-kratio)**(beta-1.) - &
            !          (kratio**(alpha-1.))*((1.-kratio)**(beta-2.))*(beta-1.)
      
            !- Beta PDF
            do k=1+kts+kb_adj-1,min(kte,kt)
                kratio=float(k+1-kb_adj)/float(kt+1)  !-kb_adj+1)
        
                zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
            enddo
            !-- special treatment below kclcl
            DO k=klcl,kts+1,-1
                zu(k)=zu(k+1)*0.5
            enddo
            !-- smooth section
            IF(do_smooth) then
                zul(kts+1)=zu(kts+1)*0.1
                ! do k=kts+2,maxloc(zu,1)
                do k=kts+2,maxloc_acc(zu)
                    zul(k)=(zu(k-1)+zu(k))*0.5
                enddo
                ! do k=kts+1,maxloc(zu,1)
                do k=kts+1,maxloc_acc(zu)
                    zu(k)=(zul(k)+zu(k))*0.5
                enddo
            ENDIF
            zu(kts)=0.
      
        !---------------------------------------------------------
        ELSEIF(draft(1:4) == "DOWN" ) then
            ! IF(trim(cumulus) == 'mid' ) beta =6.5 ! P10
            IF(cumulus(1:3) == 'mid' ) beta =6.5 ! P10

           !IF(trim(cumulus) == 'mid' ) beta =2.5 ! P9_P6
            
            ! IF(trim(cumulus) == 'deep') beta =2.5
            IF(cumulus(1:4) == 'deep') beta =2.5
           
            !IF(trim(cumulus) == 'deep') beta =6.5 ! DNDP
            !~ if(xland  > 0.90 ) then !- over water
                  !~ beta = 2.5 !ocean - lower level
            !~ else!- over land
                  !~ beta = 1.5 !land
            !~ endif
            !
            if(xland  < 0.90 ) then !- over land
                !kb_adj=kb+2
                !kb_adj=kt-2
                !kb_adj=max(kts+2,kb_adj)
                !kb_adj=min(kb_adj,kt)
                pmaxzu= 0.5 * po_cup(kt) + 0.5*psur
                ! kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)!;print*,"1=",kb_adj,po_cup(kb_adj)
                kb_adj=minloc_acc(abs(po_cup(kts:kt)-pmaxzu))
                !print*,"LAND",kb_adj,pmaxzu
            else
                kb_adj=kb+2
                !kb_adj=kt-2
                kb_adj=max(kts+2,kb_adj)
                kb_adj=min(kb_adj,kt)
                pmaxzu= 0.5 * po_cup(kt) + 0.5*psur
                ! kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)!;print*,"1=",kb_adj,po_cup(kb_adj)
                kb_adj=minloc_acc(abs(po_cup(kts:kt)-pmaxzu))
                !print*,"OCEAN",kb_adj,pmaxzu
            endif
      
            !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
            alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
            
            do k=kts+1,min(kt+1,ktf)
                kratio= float(k)/float(kt+1)
                zu(k)= kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
                !write(0,*)k,zuh(k)
            enddo
            !-- smooth section
            IF(do_smooth) then
                zul(kts+1)=zu(kts+1)*0.1
                wgty=0.
                ! do k=kts+2,maxloc(zu,1)
                do k=kts+2,maxloc_acc(zu)
                    ! wgty=wgty+1./(float(max(2,maxloc(zu,1)))-1.)
                    wgty=wgty+1./(float(max(2,maxloc_acc(zu)))-1.)
                    wgty=0.5
                    !print*,"zD1=",k,zu(k),zul(k-1),wgty,zul(k-1)*(1.-wgty)+ zu(k)*wgty
                    zul(k)=zul(k-1)*(1.-wgty)+ zu(k)*wgty
                enddo
                wgty=0.
                ! do k=kts+1,maxloc(zu,1)
                do k=kts+1,maxloc_acc(zu)
                    ! wgty=wgty+1./(float(max(2,maxloc(zu,1)))-1.)
                    wgty=wgty+1./(float(max(2,maxloc_acc(zu)))-1.)
                    wgty=0.5
                    !print*,"zD2=",k,zu(k),zul(k),wgty,zul(k)*(1.-wgty)+ zu(k)*wgty
                    zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
                enddo
            ENDIF
            zu(kts)=0.
        ENDIF
      
        ! ---------------------------------------------------------
      
        IF( maxval(zu(kts:min(kte,kt+1)),1) <= 0.0) then
            zu=0.0
            ierr=51 !ierr(i)=51
        ELSE
            !- normalize ZU
            zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ (1.e-9+maxval(zu(kts:min(kte,kt+1)),1))
        ENDIF
        
    end SUBROUTINE get_zu_zd_pdf

    SUBROUTINE get_lateral_massflux(itf,ktf, its,ite, kts,kte                             &
        ,ierr,ktop,zo_cup,zuo,cd,entr_rate_2d                 &
        ,up_massentro, up_massdetro ,up_massentr, up_massdetr &
        ,draft,kbcon,k22,up_massentru,up_massdetru,lambau)

        implicit none
        ! character *(*), intent (in) :: draft
        character *(10), intent (in) :: draft
        integer, intent(in):: itf,ktf, its,ite, kts,kte
        integer, intent(in)   , dimension(its:ite)            :: ierr,ktop,kbcon,k22
        real,    intent(in)   , dimension(its:ite), OPTIONAL  :: lambau
        real,    intent(in)   , dimension(its:ite,kts:kte) :: zo_cup,zuo
        real,    intent(inout), dimension(its:ite,kts:kte) :: cd,entr_rate_2d
        real,    intent(  out), dimension(its:ite,kts:kte) :: up_massentro, up_massdetro&
                                        ,up_massentr,  up_massdetr
        real,    intent(  out), dimension(its:ite,kts:kte),  &
                            OPTIONAL :: up_massentru,up_massdetru
        !
        !-- local vars
        integer :: i,k, incr1,incr2,turn, zuo_maxloc
        real :: dz
        logical, parameter :: SMOOTH = .FALSE.
        integer, parameter :: mass_U_option = 1  
        !---
        if(draft(1:4) == 'deep' ) then
            incr1=1
            incr2=1
        elseif(draft(1:7) == 'shallow' .or. draft(1:3) == 'mid') then
            incr1=1
            incr2=1
        else
            stop 'unknow draft in get_lateral_massflux'
        endif

!$acc kernels
        up_massentro(:,:)=0.
        up_massdetro(:,:)=0.
        if(present(up_massentru) .and. present(up_massdetru))then
            up_massentru(:,:)=0.
            up_massdetru(:,:)=0.
        endif
!$acc end kernels

!$acc kernels
        do i=its,itf
            if(ierr(i)/= 0)cycle
            zuo_maxloc = maxloc_acc(zuo(i,:))
            !-will not allow detrainment below the location of the maximum zu
            if(draft(1:7)=='shallow'.or.draft(1:3) == 'mid') cd(i,1:zuo_maxloc-2)=0.0

            !- mass entrainment and detrainment are defined on model levels

            do k=kts+1,zuo_maxloc
                !=> below location of maximum value zu -> change entrainment

                dz=zo_cup(i,k)-zo_cup(i,k-1)
                up_massdetro(i,k-1)=cd(i,k-1)*dz*zuo(i,k-1)
                up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)+up_massdetro(i,k-1)

                up_massentro(i,k-1)=max(up_massentro(i,k-1),0.0)
                !-- check dd_massdetro in case of dd_massentro has been changed above
                up_massdetro(i,k-1)=-zuo(i,k)+zuo(i,k-1)+up_massentro(i,k-1)

                if(zuo(i,k-1).gt.0.) then
                    cd(i,k-1)          =up_massdetro(i,k-1)/(dz*zuo(i,k-1))
                    entr_rate_2d(i,k-1)=up_massentro(i,k-1)/(dz*zuo(i,k-1))
                endif

            enddo

            !=================
            IF(SMOOTH .and. draft(1:4) == 'deep') THEN ! only for deep updraft
                !---smoothing the transition zone (from maxloc(zu)-1 to maxloc(zu)+1)
                do k=zuo_maxloc-1,zuo_maxloc+1
                    dz=zo_cup(i,k)-zo_cup(i,k-1)
                    up_massentro(i,k-1)=0.5*(entr_rate_2d(i,k-1)*dz*zuo(i,k-1)+up_massentro(i,k-2))
                    up_massdetro(i,k-1)=zuo(i,k-1)+up_massentro(i,k-1)-zuo(i,k)

                    if(up_massdetro(i,k-1).lt.0.)then
                        up_massdetro(i,k-1)=0.
                        up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)
                        if(zuo(i,k-1).gt.0.) entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
                    endif

                    if(zuo(i,k-1).gt.0.)cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
                enddo

                do k=zuo_maxloc,zuo_maxloc+2
                    dz=zo_cup(i,k)-zo_cup(i,k-1)

                    up_massdetro(i,k-1)=0.5*(cd(i,k-1)*dz*zuo(i,k-1)+up_massdetro(i,k-2))
                    up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)+up_massdetro(i,k-1)
                    if(up_massentro(i,k-1).lt.0.)then
                        up_massentro(i,k-1)=0.
                        up_massdetro(i,k-1)=zuo(i,k-1)-zuo(i,k)
                        if(zuo(i,k-1).gt.0.) cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
                    endif
                    if(zuo(i,k-1).gt.0.) &
                        entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
                enddo
            !-----end of the transition zone
            ENDIF
            !=================

            do k=zuo_maxloc+incr1,ktop(i)+incr2
                !=> above location of maximum value zu -> change detrainment
                dz=zo_cup(i,k)-zo_cup(i,k-1)

                up_massentro(i,k-1)=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
                !-special treatment for ktop (only for shallow)
                if(k-1==ktop(i))up_massentro(i,k-1)=0.0

                up_massdetro(i,k-1)=zuo(i,k-1)+up_massentro(i,k-1)-zuo(i,k)
                up_massdetro(i,k-1)=max(up_massdetro(i,k-1),0.0)
                !-- check up_massentro in case of dd_up_massdetro has been changed above
                up_massentro(i,k-1)=-zuo(i,k-1)+up_massdetro(i,k-1)+zuo(i,k)

                if(zuo(i,k-1).gt.0.) then
                    cd(i,k-1)          =up_massdetro(i,k-1)/(dz*zuo(i,k-1))
                    entr_rate_2d(i,k-1)=up_massentro(i,k-1)/(dz*zuo(i,k-1))
                endif
            enddo

            do k=kts,kte
                up_massentr(i,k)=up_massentro(i,k)
                up_massdetr(i,k)=up_massdetro(i,k)
            enddo
            IF(present(up_massentru) .and. present(up_massdetru))THEN
                if(mass_U_option==1) then 
                    do k=kts+1,kte
                        !--       for weaker mixing
                        up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
                        up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
                        !--       for stronger mixing
                        ! up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massentro(i,k-1)
                        ! up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massentro(i,k-1)
                    enddo
                else
                    turn=zuo_maxloc
                    do k=kts+1,turn
                        up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massentro(i,k-1)
                        up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massentro(i,k-1)
                    enddo
                    do k=turn+1,kte
                        up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
                        up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
                    enddo
                endif	 
            ENDIF    
            do k=ktop(i)+1,kte
                cd          (i,k)=0.
                entr_rate_2d(i,k)=0.
            enddo

        ENDDO
!$acc end kernels
    END SUBROUTINE get_lateral_massflux

    SUBROUTINE cup_up_moisture(name,klcl,ierr,z_cup,qc,qrc,pw,pwav,hc,tempc,xland,   &
        po,p_cup,kbcon,ktop,cd,dby,clw_all,&
        t_cup,q,GAMMA_cup,zu,qes_cup,k22,qe_cup,         &
        ZQEXEC,use_excess,ccn,rho,                       &
        up_massentr,up_massdetr,psum,psumh,c1d,          &
        itest,itf,ktf,ipr,jpr,its,ite, kts,kte           )

        implicit none
        real, parameter :: bdispm = 0.366       !berry--size dispersion (maritime)
        real, parameter :: bdispc = 0.146       !berry--size dispersion (continental)
        real, parameter :: T_BF  = 268.16   &
                    , T_ice = 250.16
        !
        !  on input
        !

        ! only local  dimensions are need as of now in this routine

        integer                                                           &
        ,intent (in   )                   ::                           &
                                    use_excess,itest,itf,ktf,            &
                                    its,ite,ipr,jpr, kts,kte
        ! cd= detrainment function
        ! q = environmental q on model levels
        ! qe_cup = environmental q on model cloud levels
        ! qes_cup = saturation q on model cloud levels
        ! dby = buoancy term
        ! cd= detrainment function
        ! zu = normalized updraft mass flux
        ! gamma_cup = gamma on model cloud levels
        !
        real,    dimension (its:ite,kts:kte)                              &
            ,intent (in   )                   ::                           &
            t_cup,po,p_cup,rho,q,zu,gamma_cup,qe_cup,hc,                   &
            up_massentr,up_massdetr,dby,qes_cup,z_cup,cd,c1d
        real,    dimension (its:ite)                                      &
            ,intent (in   )                   ::                           &
            zqexec,xland
        integer, dimension (its:ite)                                      &
            ,intent (in   )                   ::                           &
            kbcon,ktop,k22,klcl
        !
        ! input and output
        !

        ! ierr error value, maybe modified in this routine

        integer, dimension (its:ite)                                      &
            ,intent (inout)                   ::                           &
            ierr
        ! character *(*), intent (in)        ::                           &
        character *(10), intent (in)        ::                           &
            name
        ! qc = cloud q (including liquid water) after entrainment
        ! qrch = saturation q in cloud
        ! qrc = liquid water content in cloud after rainout
        ! pw = condensate that will fall out at that level
        ! pwav = totan normalized integrated condensate (I1)
        ! c0 = conversion rate (cloud to rain)

        real,    dimension (its:ite,kts:kte)                              &
            ,intent (out  )                   ::                           &
            qc,qrc,pw,clw_all,tempc
        real,    dimension (its:ite,kts:kte) ::                           &
            qch,qrcb,pwh,clw_allh
        ! real,    dimension (its:ite)         ::                           &
        !     pwavh
        real,    dimension (its:ite)                                      &
            ,intent (out  )                   ::                           &
            pwav,psum,psumh
        real,    dimension (its:ite)                                      &
            ,intent (in  )                   ::                           &
            ccn
        !
        !  local variables in this routine
        !
        integer                              ::                           &
            iounit,iprop,i,k,k1,k2,n,nsteps,rk
        real                                 ::                           &
            dp,rhoc,dh,qrch,c0,dz,radius,berryc0,q1,berryc
        real :: qaver,denom,aux,cx0,qrci,step,cbf,qrc_crit_BF,min_liq,xexp,qavail
        integer, dimension (its:ite)   ::     start_level
        real delt

!$acc data create(qch,qrcb,pwh,clw_allh, start_level) copyin(name)


!$acc kernels
        start_level =0
!$acc end kernels
        !
        !--- no precip for small clouds
        if(name(1:7).eq.'shallow')then
            c0=0.0
!$acc kernels
            start_level(:)=klcl(:)!kbcon
!$acc end kernels
        endif
        if(name(1:3).eq.'mid' )then
            c0=c0_mid
!$acc kernels
            start_level(:)=klcl(:)!kbcon(:)
!$acc end kernels
        endif
        if(name(1:4).eq.'deep' )then
            c0=c0_deep
!$acc kernels
            start_level(:)=klcl(:)!k22(:)
!$acc end kernels
        endif

!$acc parallel loop vector present(pwav, psum, psumh)
        do i=its,itf
            pwav (i)=0.
            ! pwavh(i)=0.
            psum (i)=0.
            psumh(i)=0.
        enddo
!$acc end parallel

!$acc parallel 
!$acc loop gang vector collapse(2)
        do k=kts,kte
            do i=its,itf
                pw      (i,k)=0.
                pwh     (i,k)=0.
                qc      (i,k)=0.
                qrc     (i,k)=0.
                qrcb    (i,k)=0.
                clw_all (i,k)=0.
                clw_allh(i,k)=0.
                tempc   (i,k)=t_cup(i,k)
            enddo
        enddo
!$acc loop vector
        do i=its,itf
            if(ierr(i) /= 0)cycle
            do k=kts,ktf
                qc (i,k)=qe_cup(i,k)
                qch(i,k)=qe_cup(i,k)
            enddo
        enddo
!$acc end parallel


!$acc kernels
        do i=its,itf
            if(ierr(i).eq.0)then
                call get_cloud_bc(name,kts,kte,ktf,xland(i),po(i,kts:kte),qe_cup (i,kts:kte),qaver,k22(i))
                qaver = qaver + float(use_excess)*zqexec(i)

                do k=kts,start_level(i)!k22(i)
                    qc (i,k)=qaver
                    qch(i,k)=qc(i,k)
                
                enddo
            endif
        enddo
!$acc end kernels

!$acc kernels
        DO 100 i=its,itf
            IF (ierr(i)  /= 0) cycle

            DO k=start_level(i)+1,ktop(i) + 1     !k22(i)+1,ktop(i)

                DZ=Z_cup(i,K)-Z_cup(i,K-1)
                !
                !--- saturation  in cloud, this is what is allowed to be in it
                !
                QRCH = qes_cup(I,K)+(1./XLV)*(GAMMA_cup(i,k) &
                    /(1.+GAMMA_cup(i,k)))*DBY(I,K)

                !-    1. steady state plume equation, for what could
                !-       be in cloud without condensation
                denom =  (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
                if(denom > 0.) then
                    qc (i,k)=  ( qc (i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1) +   &
                                                        up_massentr(i,k-1)* q (i,k-1)     &
                                )/ denom
                else
                    qc (i,k)=    qc (i,k-1)
                endif
                !-- updraft temp
                tempc(i,k) = (1./cp)*(hc(i,k)-g*z_cup(i,k)-xlv*qc(i,k))

                !------- Total condensed water before rainout
                clw_all(i,k)=max(0.,QC(I,K)-QRCH)

                if(name(1:7).eq.'shallow')then
                    qrc(i,k) = clw_all(i,k)
                    qc (i,k) = qrc(i,k)+min(qc(i,k),qrch)
                    pwav(i)=0.
                    psum(i)=psum(i)+clw_all(i,k)*zu(i,k) *dz
                    cycle
                ENDIF

                IF (autoconv == 1 ) then
                    !--- ok  exp a1_5
                    cx0     = (c1d(i,k)+c0)*DZ!*zu(i,k)
                    qrc(i,k)= clw_all(i,k)/(1.+cx0)
                    PW (i,k)= cx0*max(0.,QRC(I,K) -qrc_crit)! units kg[rain]/kg[air]
                    !--- convert PW to normalized PW
                    PW (i,k)=PW(i,k)*zu(i,k)

                    !--- test exp a1_6 (not ok)
                    !.. cx0 = (c1d(i,k)+c0)*DZ!*zu(i,k)
                    !.. qavail = max(0.,clw_all(i,k)-qrc_crit)
                    !.. !- liq water after rainout
                    !.. qavail =  qavail/(1+cx0)
                    !.. pw (i,k)= cx0*qavail
                    !.. QRC(i,k)= clw_all(i,k)-pw(i,k)
                    !.. !--- convert PW to normalized PW
                    !.. PW (i,k)=PW(i,k)*zu(i,k)

                ELSEIF (autoconv == 11 ) then
                    IF(clw_all(i,k) <= qrc_crit) then
                        QRC(I,K)= clw_all(i,k)
                        PW(i,k) = 0.
                    ELSE
                        cx0     = (c1d(i,k)+c0)*DZ!*zu(i,k)
                        QRC(i,k)= (clw_all(i,k))/(1.+cx0)
                        PW (i,k)= cx0*(QRC(i,k))!*zu(i,k)! units kg[rain]/kg[air]
                        !--- convert PW to normalized PW
                        PW(i,k)=PW(i,k)*zu(i,k)
                    ENDIF
                ELSEIF (autoconv == 3 ) then
                    !c0=.001
                    DELT=-20.
                    if(t_cup(i,k) > 273.16 + DELT) then
                        aux = 1.
                    else
                        aux = 1. * exp(0.03* (t_cup(i,k) - 273.16 + DELT))
                    endif
  
                    cx0 = aux*(c1d(i,k)+c0)*DZ!*zu(i,k) !exp B1a3x2
                    !-v1 -  do not conserve total water:     QRCH+PW (i,k)+QRC(i,k)
                    !.. QRC(i,k)=max(0.,  QC(i,k)- QRCH + cx0*DZ*qrc_crit)/(1.+cx0*DZ)
                    !.. PW (i,k)=max(0., QRC(i,k) -qrc_crit)*cx0*DZ
                    !.. !--- convert PW to normalized PW
                    !.. PW (i,k) = PW(i,k)*zu(i,k)
                    !-v2 - conserves total water:     QRCH+PW (i,k)+QRC(i,k)
                    QRC(i,k)=(clw_all(i,k))/(1.+cx0)
                    PW (i,k)=max(0., QRC(i,k) -qrc_crit)*cx0
                    !print*,"cons=",k,PW (i,k)+QRCH+QRC(i,k),qc(i,k)
                    PW (i,k) = PW(i,k)*zu(i,k)
                ELSEIF (autoconv == 4 ) then

                    qrc(i,k) = clw_all(i,k)!- initial cloud water
                    min_liq  = (xland(i)*0.3+(1.-xland(i))*0.5)*1.e-3

                    if(qrc(i,k)>min_liq) then
                        cbf = 1.
                        if(tempc(i,k) < T_BF) cbf=1.+0.5*sqrt (min(T_BF-tempc(i,k),T_BF-T_ICE))
                        !cx0 = 0.002*cbf
                        cx0 = (c1d(i,k)+c0)*cbf !use this: c1d(i,k)+c0
                        !print*,'cbf=',k,cx0,cbf,dz
                        qrc_crit_BF=qrc_crit/cbf
                        rk = 3 ! or 2
                        xexp = 2.
                        IF(clw_all(i,k)> 1.e-10) then
                            step = cx0*dz!*zu(i,k)
                            DO n=rk,1,-1
                                aux     = qrc(i,k)/qrc_crit_BF; aux=max(aux,2.5)
                                pw (i,k)= auto_rk(n,step,aux,xexp,qrc(i,k))
                                qrc(i,k)= max(clw_all(i,k) - pw(i,k), min_liq)
                                !if(n==1)print*,"rk3N=",k,qrc(i,k),PW(i,k),QC(I,K)-QRCH
                            ENDDO
                            !- for checking the solution
                            !.. aux  = qrc(i,k)/qrc_crit_BF; aux=max(aux,2.5)
                                    !.. pw (i,k)=auto_rk(1,step,aux,xexp,qrc(i,k))
                            !.. print*,"rk3F=",k,dz,qrc(i,k),PW(i,k)
                                    !--- convert PW to normalized PW
                            PW (i,k) = PW(i,k)*zu(i,k)
                        ENDIF
                    ENDIF
                ENDIF
                !- total water (vapor + condensed) in updraft after the rainout
                !print*,"qc=",k,qc(i,k)*1000,(QRC(I,K)+min(QC(i,k),QRCH)   )*1000.,QRC(I,K)  *1000.,pw(i,k)*1000.
                QC(I,K)=QRC(I,K)+min(QC(i,k),QRCH)

                !--- integrated normalized condensates
                PWAV(I)=PWAV(I)+PW(I,K)
                psum(I)=psum(I)+clw_all(I,K)*zu(i,k) *dz

            ENDDO
            !.. ELSEIF( autoconv == 2) then
            !.. STOP " BERRY NOT READY YET"

            !.. ENDIF
100     CONTINUE
!$acc end kernels
        !- get back water vapor qc
!$acc parallel loop gang
        do i=its,itf
            if(ierr(i).eq.0)then
!$acc loop vector
                do k=kts,ktop(i)+1
                    qc(i,k)= qc(i,k)-qrc(i,k)
                    !if(qc(i,k) < 0.)stop " qc negative"
                enddo
            endif
        enddo
!$acc end parallel

!$acc end data

    END SUBROUTINE cup_up_moisture

    subroutine cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd ,z,z_cup,zu,dby,GAMMA_CUP,t_cup &
        ,tempco,qco,qrco,qo,klcl,kbcon,ktop &
        ,ierr,itf,ktf,its,ite, kts,kte)

        implicit none
        real, parameter :: ctea=1./3. ,cteb=2., visc=2000., eps=0.622
        integer,intent (in   )              ::  itf,ktf,its,ite, kts,kte
        real,    dimension (its:ite,kts:kte) ,intent (in   )  ::  &
            z,z_cup,zu,gamma_cup,t_cup,dby,entr_rate_2d,cd,tempco,qco,qrco,qo

        integer, dimension (its:ite)         ,intent (in   )  ::  &
            klcl,kbcon,ktop
        real,    dimension (its:ite)         ,intent (in   )  ::  &
            zws

        ! input and output
        integer, dimension (its:ite)        ,intent (inout) :: ierr
        real,    dimension (its:ite,kts:kte),intent (out  ) ::  vvel2d
        real,    dimension (its:ite        ),intent (out  ) ::  vvel1d
        !
        !  local variables in this routine
        integer                             ::  i,k,k1,nvs
        real                                ::  dz,BU,dw2,dw1,kx,dz1m,Tv,Tve,vs
        real   , parameter :: f=2., C_d=0.506, gam=0.5, beta=1.875
        logical, parameter :: smooth=.true.
        integer, parameter :: n_smooth=3

!$acc kernels
        !-- initialize arrays to zero.
        vvel1d = 0.0
        vvel2d = 0.0
!$acc end kernels

!$acc parallel loop gang
        do i=its,itf
            ! !-- initialize arrays to zero.
            ! vvel1d(i  ) = 0.0
            ! vvel2d(i,:) = 0.0

            if(ierr(i) /= 0) cycle
            vvel2d(i,kts:kbcon(i))= max(1.,zws(i)**2)

!$acc loop seq
            loop0:  do k= kbcon(i),ktop(i)

                dz=z_cup(i,k+1)-z_cup(i,k)

                Tve= 0.5* ( t_cup (i,k  )*(1.+(qo (i,k  )/eps)/(1.+qo (i,k  ))) + &
                    t_cup (i,k+1)*(1.+(qo (i,k+1)/eps)/(1.+qo (i,k+1))))

                Tv = 0.5* ( tempco(i,k  )*(1.+(qco(i,k  )/eps)/(1.+qco(i,k  ))) + &
                    tempco(i,k+1)*(1.+(qco(i,k+1)/eps)/(1.+qco(i,k+1)) ))

                !--- orig
                !BU = g*( (Tv-Tve)/Tve - 0.50*(qrco(i,k+1)+qrco(i,k) ))
                BU = g*( (Tv-Tve)/Tve - 0.25*(qrco(i,k+1)+qrco(i,k) ))

                dw1 = 2./(f*(1.+gam)) * BU * dz

                !--- orig
                !kx  = (1.+beta*C_d)*max(entr_rate_2d(i,k),cd(i,k))*dz
                kx  =               max(entr_rate_2d(i,k),cd(i,k))*dz
                dw2 =  (vvel2d(i,k)) -2.*kx * (vvel2d(i,k))

                vvel2d(i,k+1)=(dw1+dw2)/(1.+kx)

                if( vvel2d(i,k+1)< 0.) then
                    vvel2d(i,k+1) = 0.5* vvel2d(i,k)
                endif

            enddo loop0
        enddo
!$acc end parallel


        if(smooth) then
!$acc parallel loop gang
            do i=its,itf
                if(ierr(i) /= 0)cycle
!$acc loop seq
                do k=kts,ktop(i)-2
                    nvs=0; vs =0.
                    do k1 = max(k-n_smooth,kts),min(k+n_smooth,ktf)
                        nvs = nvs + 1
                        vs =  vs + vvel2d(i,k1)
                    enddo
                    vvel2d(i,k) = vs/(1.e-16+float(nvs))
                enddo 
            enddo
!$acc end parallel 
        endif

        !-- convert to vertical velocity
!$acc kernels
        do i=its,itf
            if(ierr(i) /= 0)cycle
            vvel2d(i,:)= sqrt(vvel2d(i,:))

            !-- sanity check
            where(vvel2d(i,:) < 1. ) vvel2d(i,:)=1.
            where(vvel2d(i,:) > 20.) vvel2d(i,:)=20.

            !-- get the column average vert velocity
            do k= kbcon(i),ktop(i)
                dz=z_cup(i,k+1)-z_cup(i,k)
                vvel1d(i)=vvel1d(i)+vvel2d(i,k)*dz
            enddo
            vvel1d(i)=vvel1d(i)/(z_cup(i,ktop(i)+1)-z_cup(i,kbcon(i)))
            vvel1d(i)=max(1.,vvel1d(i))
        enddo
!$acc end kernels

    end subroutine cup_up_vvel

    SUBROUTINE cup_up_cape(aa0,z,zu,dby,GAMMA_CUP,t_cup,  &
        k22,kbcon,ktop,ierr,tempco,qco,qrco, qo_cup,   &
        itf,ktf,its,ite, kts,kte      )

        IMPLICIT NONE
        integer ,intent (in   )                   ::        &
            itf,ktf, its,ite, kts,kte

        ! aa0 = dummy array for CAPE (total cape)
        ! gamma_cup = gamma on model cloud levels
        ! t_cup = temperature (Kelvin) on model cloud levels
        ! dby = buoancy term
        ! zu= normalized updraft mass flux
        ! z = heights of model levels
        ! ierr = error value, maybe modified in this routine
        ! tempco = in-cloud temperature (Kelvin) on model cloud levels
        ! qco    = in-cloud water vapor mixing ratio on model cloud levels
        ! qo_cup = environ water vapor mixing ratio on model cloud levels
        ! qrco   = in-cloud liquid water mixing ratio on model cloud levels

        real,    dimension (its:ite,kts:kte) ,intent (in   )    ::       &
            z,zu,gamma_cup,t_cup,dby,tempco,qco,qrco, qo_cup
        integer, dimension (its:ite)         ,intent (in   )    ::       &
            k22,kbcon,ktop
        !
        ! input and output
        integer, dimension (its:ite)         ,intent (inout)    ::       &
            ierr
        real,    dimension (its:ite)         ,intent (out  )    ::       &
            aa0
        !
        !  local variables in this routine
        integer       ::   i,k
        real          ::   dz,daa0
        !

!$acc parallel
        aa0(:)=0.
!$acc loop gang   
        DO i=its,itf
            IF(ierr(i) == 0) then
!$acc loop vector private(dz, daa0)
                DO k=max(k22(i),kts+1),ktop(i)
                    dz=z(i,k)-z(i,k-1)
                    daa0=g*dz*(   (tempco(i,k)*(1.+0.608*qco   (i,k))  - t_cup(i,k)*(1.+0.608*qo_cup(i,k)))&
                                / (t_cup (i,k)*(1.+0.608*qo_cup(i,k))) &
                            )
!$acc atomic update
                    aa0(i)=aa0(i)+daa0
!$acc end atomic
                    !~ print*,"cape",k,AA0(I),tempco(i,k),t_cup(i,k), qrco  (i,k)
                ENDDO
            ENDIF
        ENDDO
!$acc end parallel
    END SUBROUTINE cup_up_cape

    SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
        us, vs,u_cup,v_cup,                                   &
        hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,tsur,        &
        ierr,z1,itf,ktf,its,ite, kts,kte                      &
        ,host_model )

        implicit none
        character(len=*),intent(in),optional :: host_model
        integer ,intent (in   )              :: itf,ktf, its,ite, kts,kte
        !
        ! ierr error value, maybe modified in this routine
        ! q           = environmental mixing ratio
        ! q_cup       = environmental mixing ratio on cloud levels
        ! qes         = environmental saturation mixing ratio
        ! qes_cup     = environmental saturation mixing ratio on cloud levels
        ! t           = environmental temp
        ! t_cup       = environmental temp on cloud levels
        ! p           = environmental pressure
        ! p_cup       = environmental pressure on cloud levels
        ! z           = environmental heights
        ! z_cup       = environmental heights on cloud levels
        ! he          = environmental moist static energy
        ! he_cup      = environmental moist static energy on cloud levels
        ! hes         = environmental saturation moist static energy
        ! hes_cup     = environmental saturation moist static energy on cloud levels
        ! gamma_cup   = gamma on cloud levels
        ! psur        = surface pressure
        ! z1          = terrain elevation
        !
        !
        real,    dimension (its:ite,kts:kte)                              &
            ,intent (in   )                   ::                           &
            qes,q,he,hes,z,p,t,us, vs
        real,    dimension (its:ite,kts:kte)                              &
            ,intent (out  )                   ::                           &
            qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,u_cup,v_cup
        real,    dimension (its:ite)                                      &
            ,intent (in   )                   ::                           &
            psur,z1,tsur
        integer, dimension (its:ite)                                      &
            ,intent (inout)                   ::                           &
            ierr
        !
        !  local variables in this routine
        !
        integer                              ::      i,k
        real :: p1,p2,c1,c2

!$acc kernels present(qes_cup, q_cup, hes_cup, he_cup, &
!$acc                 z_cup, p_cup, t_cup, gamma_cup, &
!$acc                 u_cup, v_cup)
        qes_cup  =0.
        q_cup    =0.
        hes_cup  =0.
        he_cup   =0.
        z_cup    =0.
        p_cup    =0.
        t_cup    =0.
        gamma_cup=0.
        u_cup    =0.
        v_cup    =0.
!$acc end kernels



        !      IF(.NOT. present(host_model)) THEN
        IF( present(host_model) .and. trim(host_model)=='OLD_GEOS5') THEN
!$acc parallel present(qes,q,he,hes,z,p,t,us,vs,qes_cup,q_cup,he_cup,&
!$acc                  hes_cup,z_cup,p_cup,gamma_cup,t_cup,u_cup,v_cup, &
!$acc                  psur,z1,tsur,ierr)

!$acc loop gang vector collapse(2)
            do k=kts+1,ktf
                do i=its,itf
                    if(ierr(i).eq.0)then
                        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
                        q_cup  (i,k)=.5*(q  (i,k-1)+q  (i,k))
                        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
                        he_cup (i,k)=.5*(he (i,k-1)+he (i,k))
                        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
                        z_cup  (i,k)=.5*(z(i,k-1)+z(i,k))
                        p_cup  (i,k)=.5*(p(i,k-1)+p(i,k))
                        t_cup  (i,k)=.5*(t(i,k-1)+t(i,k))
                        gamma_cup(i,k)=(xlv/cp)*(xlv/(rv*t_cup(i,k) &
                                        *t_cup(i,k)))*qes_cup(i,k)
                        u_cup  (i,k)=.5*(us(i,k-1)+us(i,k))
                        v_cup  (i,k)=.5*(vs(i,k-1)+vs(i,k))                            
                    endif
                enddo
            enddo

!$acc loop vector
            do i=its,itf
                if(ierr(i).eq.0)then
                    qes_cup(i,1)=qes(i,1)
                    q_cup(i,1)=q(i,1)
                    !       hes_cup(i,1)=hes(i,1)
                    !       he_cup(i,1)=he(i,1)
                    hes_cup(i,1)=g*z1(i)+1004.*t(i,1)+2.5e6*qes(i,1)
                    he_cup (i,1)=g*z1(i)+1004.*t(i,1)+2.5e6*q  (i,1)
                    !       z_cup(i,1)=.5*(z(i,1)+z1(i))
                    !       p_cup(i,1)=.5*(p(i,1)+psur(i))
                    z_cup(i,1)=z1(i)
                    p_cup(i,1)=psur(i)
                    t_cup(i,1)=t(i,1)
                    gamma_cup(i,1)=xlv/cp*(xlv/(rv*t_cup(i,1) &
                                    *t_cup(i,1)))*qes_cup(i,1)
                    u_cup(i,1)=us(i,1)
                    v_cup(i,1)=vs(i,1)
                
                endif
            enddo
!$acc end parallel

        ELSE
            IF(trim(host_model)=='XXXGEOS5') THEN
                print*,'Branch gets exercised 1'
!$acc parallel present(qes,q,he,hes,z,p,t,us,vs,qes_cup,q_cup,he_cup,&
!$acc                  hes_cup,z_cup,p_cup,gamma_cup,t_cup,u_cup,v_cup, &
!$acc                  psur,z1,tsur,ierr)

!$acc loop vector
                do i = its,itf
                    if(ierr(i) > 0) cycle
                    
                    qes_cup(i,1)=qes(i,1)
                    q_cup  (i,1)=q(i,1)
                    hes_cup(i,1)=g*z1(i)+1004.*t(i,1)+2.5e6*qes(i,1)
                    he_cup (i,1)=g*z1(i)+1004.*t(i,1)+2.5e6*q  (i,1)

                    z_cup  (i,1)=z1(i)
                    p_cup  (i,1)=psur(i)
                    t_cup  (i,1)=t(i,1)
                    u_cup  (i,1)=us(i,1)
                    v_cup  (i,1)=vs(i,1)
                    gamma_cup(i,1)=xlv/cp*(xlv/(rv*t_cup(i,1)*t_cup(i,1)))*qes_cup(i,1)
                enddo

!$acc loop seq
                do k = kts,ktf-1
!$acc loop vector
                    do i = its,itf
                        if(ierr(i) > 0) cycle
                        z_cup  (i,k+1)=2.0*z  (i,k)-z_cup  (i,k)
                        p_cup  (i,k+1)=2.0*p  (i,k)-p_cup  (i,k)
                        t_cup  (i,k+1)=2.0*t  (i,k)-t_cup  (i,k)
                        u_cup  (i,k+1)=2.0*us (i,k)-u_cup  (i,k)
                        v_cup  (i,k+1)=2.0*vs (i,k)-v_cup  (i,k)
                        q_cup  (i,k+1)=2.0*q  (i,k)-q_cup  (i,k)
                        he_cup (i,k+1)=2.0*he (i,k)-he_cup (i,k)

                        qes_cup(i,k+1)=2.0*qes(i,k)-qes_cup(i,k)
                        hes_cup(i,k+1)=2.0*hes(i,k)-hes_cup(i,k)
    
                        if(he_cup(i,k+1).gt.hes_cup(i,k+1))he_cup(i,k+1)=hes_cup(i,k+1)

                        gamma_cup(i,k+1)=(xlv/cp)*(xlv/(rv*t_cup(i,k+1) &
                                    *t_cup(i,k+1)))*qes_cup(i,k+1)
                    enddo
                enddo

!$acc end parallel

                ! do i=its,itf
                !     if(ierr(i) > 0) cycle
                    
                !     qes_cup(i,1)=qes(i,1)
                !     q_cup  (i,1)=q(i,1)
                !     hes_cup(i,1)=g*z1(i)+1004.*t(i,1)+2.5e6*qes(i,1)
                !     he_cup (i,1)=g*z1(i)+1004.*t(i,1)+2.5e6*q  (i,1)

                !     z_cup  (i,1)=z1(i)
                !     p_cup  (i,1)=psur(i)
                !     t_cup  (i,1)=t(i,1)
                !     u_cup  (i,1)=us(i,1)
                !     v_cup  (i,1)=vs(i,1)
                !     gamma_cup(i,1)=xlv/cp*(xlv/(rv*t_cup(i,1)*t_cup(i,1)))*qes_cup(i,1)
    
                !     do k=kts,ktf-1
                !         z_cup  (i,k+1)=2.0*z  (i,k)-z_cup  (i,k)
                !         p_cup  (i,k+1)=2.0*p  (i,k)-p_cup  (i,k)
                !         t_cup  (i,k+1)=2.0*t  (i,k)-t_cup  (i,k)
                !         u_cup  (i,k+1)=2.0*us (i,k)-u_cup  (i,k)
                !         v_cup  (i,k+1)=2.0*vs (i,k)-v_cup  (i,k)
                !         q_cup  (i,k+1)=2.0*q  (i,k)-q_cup  (i,k)
                !         he_cup (i,k+1)=2.0*he (i,k)-he_cup (i,k)

                !         qes_cup(i,k+1)=2.0*qes(i,k)-qes_cup(i,k)
                !         hes_cup(i,k+1)=2.0*hes(i,k)-hes_cup(i,k)
    
                !         if(he_cup(i,k+1).gt.hes_cup(i,k+1))he_cup(i,k+1)=hes_cup(i,k+1)

                !         gamma_cup(i,k+1)=(xlv/cp)*(xlv/(rv*t_cup(i,k+1) &
                !                     *t_cup(i,k+1)))*qes_cup(i,k+1)
                !     enddo
                ! ENDDO
            ELSEIF(trim(host_model)=='NEW_GEOS5') THEN
!$acc parallel present(qes,q,he,hes,z,p,t,us,vs,qes_cup,q_cup,he_cup,&
!$acc                  hes_cup,z_cup,p_cup,gamma_cup,t_cup,u_cup,v_cup, &
!$acc                  psur,z1,tsur,ierr)

!$acc loop vector
                do i = its,itf
                    if(ierr(i) > 0) cycle
                    p_cup  (i,1)=psur(i)
                    z_cup  (i,1)=z1(i)
                enddo

!$acc loop seq
                do k = kts,ktf-1
!$acc loop vector
                    do i = its, itf
                        if(ierr(i) > 0) cycle
                        p_cup (i,k+1) = 2.0*p(i,k) - p_cup(i,k)
                        z_cup (i,k+1) = 2.0*z(i,k) - z_cup(i,k)
                    enddo
                enddo

!$acc loop gang vector collapse(2) private(p1,p2)
                do k = kts,ktf-1
                    do i = its,itf
                        if(ierr(i) > 0) cycle
                        p1=abs((p    (i,k+1) - p_cup(i,k+1))/(p(i,k+1)-p(i,k)))
                        p2=abs((p_cup(i,k+1) - p    (i,k  ))/(p(i,k+1)-p(i,k)))
                                                                    
                        t_cup  (i,k+1) = p1*t  (i,k) + p2*t  (i,k+1)	  	      
                                                            
                        u_cup  (i,k+1) = p1*us (i,k) + p2*us (i,k+1)
                        v_cup  (i,k+1) = p1*vs (i,k) + p2*vs (i,k+1)
                        q_cup  (i,k+1) = p1*q  (i,k) + p2*q  (i,k+1)
                        he_cup (i,k+1) = p1*he (i,k) + p2*he (i,k+1)
                                                                    
                        qes_cup(i,k+1) = p1*qes(i,k) + p2*qes(i,k+1)
                        hes_cup(i,k+1) = p1*hes(i,k) + p2*hes(i,k+1)
                                                            
                        if(he_cup(i,k+1).gt.hes_cup(i,k+1))he_cup(i,k+1)=hes_cup(i,k+1)
                                                            
                        gamma_cup(i,k+1)=(xlv/cp)*(xlv/(rv*t_cup(i,k+1) &		  	      
                                        *t_cup(i,k+1)))*qes_cup(i,k+1)
                    enddo
                enddo

                k = kts

!$acc loop vector private(p1, p2, c1, c2)
                do i = its,itf
                    p1=abs(p    (i,k  )-p_cup(i,k))
                    p2=abs(p_cup(i,k+1)-p_cup(i,k))
                                                                    
                    c1=(p1+p2)/p2
                    c2=p1/p2
                                                  
                    t_cup  (i,k) = c1*t  (i,k) - c2*t_cup(i,k+1)
                    q_cup  (i,k) = c1*q  (i,k) - c2*q_cup(i,k+1)
                                                                
                    u_cup  (i,k) = c1*us (i,1) - c2*u_cup(i,k+1)
                    v_cup  (i,k) = c1*vs (i,1) - c2*v_cup(i,k+1)
                    qes_cup(i,k) = c1*qes(i,k) - c2*qes_cup(i,k+1)
                                                                
                    hes_cup(i,k)=g*z_cup(i,k)+1004.*t_cup(i,k)+2.5e6*qes_cup(i,k)
                    he_cup (i,k)=g*z_cup(i,k)+1004.*t_cup(i,k)+2.5e6*q_cup  (i,k)
      
                    if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
      
                    gamma_cup(i,k)=xlv/cp*(xlv/(rv*t_cup(i,k)*t_cup(i,k)))*qes_cup(i,k)
                enddo

!$acc end parallel

!!$acc update host(qes_cup, q_cup, hes_cup, he_cup, z_cup, p_cup, t_cup, &
!!$acc             gamma_cup, u_cup, v_cup)

                ! do i=its,itf
                !     if(ierr(i) > 0) cycle
                !     p_cup  (i,1)=psur(i)
                !     z_cup  (i,1)=z1(i)
                !     do k=kts,ktf-1
                !         p_cup (i,k+1) = 2.0*p(i,k) - p_cup(i,k)
                !         z_cup (i,k+1) = 2.0*z(i,k) - z_cup(i,k)
                !     enddo
  
                !     ! ----------- p,T	       k+1
                !     !p1 
                !     ! ----------- p_cup,T_cup  k+1
                !     !p2 
                !     ! ----------- p,T	       k
                !     ! 
                !     ! ----------- p_cup,T_cup  k
  
                !     do k=kts,ktf-1			  	      
                !         p1=abs((p    (i,k+1) - p_cup(i,k+1))/(p(i,k+1)-p(i,k)))
                !         p2=abs((p_cup(i,k+1) - p    (i,k  ))/(p(i,k+1)-p(i,k)))
                                                                    
                !         t_cup  (i,k+1) = p1*t  (i,k) + p2*t  (i,k+1)	  	      
                                                            
                !         u_cup  (i,k+1) = p1*us (i,k) + p2*us (i,k+1)
                !         v_cup  (i,k+1) = p1*vs (i,k) + p2*vs (i,k+1)
                !         q_cup  (i,k+1) = p1*q  (i,k) + p2*q  (i,k+1)
                !         he_cup (i,k+1) = p1*he (i,k) + p2*he (i,k+1)
                                                                    
                !         qes_cup(i,k+1) = p1*qes(i,k) + p2*qes(i,k+1)
                !         hes_cup(i,k+1) = p1*hes(i,k) + p2*hes(i,k+1)
                                                            
                !         if(he_cup(i,k+1).gt.hes_cup(i,k+1))he_cup(i,k+1)=hes_cup(i,k+1)
                                                            
                !         gamma_cup(i,k+1)=(xlv/cp)*(xlv/(rv*t_cup(i,k+1) &		  	      
                !                         *t_cup(i,k+1)))*qes_cup(i,k+1)
                                                            
                !     enddo
                !     !--- surface level from X(kts) and X_cup(kts+1) determine X_cup(kts)	      
                !     k=kts
                !     p1=abs(p    (i,k  )-p_cup(i,k))
                !     p2=abs(p_cup(i,k+1)-p_cup(i,k))
                                                                    
                !     c1=(p1+p2)/p2
                !     c2=p1/p2
                                                  
                !     t_cup  (i,k) = c1*t  (i,k) - c2*t_cup(i,k+1)
                !     q_cup  (i,k) = c1*q  (i,k) - c2*q_cup(i,k+1)
                                                                
                !     u_cup  (i,k) = c1*us (i,1) - c2*u_cup(i,k+1)
                !     v_cup  (i,k) = c1*vs (i,1) - c2*v_cup(i,k+1)
                !     qes_cup(i,k) = c1*qes(i,k) - c2*qes_cup(i,k+1)
                                                                
                !     hes_cup(i,k)=g*z_cup(i,k)+1004.*t_cup(i,k)+2.5e6*qes_cup(i,k)
                !     he_cup (i,k)=g*z_cup(i,k)+1004.*t_cup(i,k)+2.5e6*q_cup  (i,k)
      
                !     if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
      
                !     gamma_cup(i,k)=xlv/cp*(xlv/(rv*t_cup(i,k)*t_cup(i,k)))*qes_cup(i,k)
                ! enddo
            ELSE

                STOP "cup_env_clev" 
            ENDIF
            
        ENDIF

!!$acc update device(qes_cup, q_cup, hes_cup, he_cup, z_cup, p_cup, t_cup, &
!!$acc             gamma_cup, u_cup, v_cup)
!!$acc end data
        !------------------------------------------------
        RETURN
        ! IF( MAPL_AM_I_ROOT()) then
        !     DO i=its,itf
        !         IF(ierr(i) == 0) then 
        !             do k=kts,kte
        !                 write(23,101) i,k,z_cup(i,k),p_cup(i,k),t_cup(i,k),u_cup(i,k),v_cup(i,k),q_cup(i,k),he_cup(i,k)
        !                 write(25,101) i,k,z    (i,k),p    (i,k),t    (i,k),us   (i,k),vs   (i,k),q    (i,k),he    (i,k)

        !                 101 format(2i3,7F15.4)
        !             ENDDO
        !         GOTO 400
        !         ENDIF
        !     ENDDO
        ! ENDIF
        ! 400 continue

    END SUBROUTINE cup_env_clev

    SUBROUTINE cup_env_clev_chem(mtp,se_chem,se_cup_chem,ierr,itf,ktf,its,ite, kts,kte)

        IMPLICIT NONE
        !-inputs
        integer ,intent (in   )                   :: itf,ktf, its,ite, kts,kte
        integer ,intent (in   )                   :: mtp
        integer, dimension (its:ite),intent (in)  :: ierr
        real, dimension (mtp,its:ite,kts:kte),intent (in)  ::   se_chem
        !-outputs
        real, dimension (mtp,its:ite,kts:kte),intent (out) ::   se_cup_chem
        !-locals
        integer ::  i,k,m
        integer, parameter ::  clev_option = 2 !- use option 2
   
        !
        IF(clev_option == 1) THEN
!$acc kernels
            !-- original version
            do i=its,itf
                if(ierr(i) /= 0) cycle
                do k=kts+1,ktf
                    se_cup_chem(1:mtp,i,k)=0.5*(se_chem(1:mtp,i,k-1)+se_chem(1:mtp,i,k))
                enddo
                se_cup_chem(1:mtp,i,kts)=se_chem(1:mtp,i,kts)
                se_cup_chem(1:mtp,i,kte)=se_chem(1:mtp,i,ktf)
            enddo
!$acc end kernels
        ELSE
            !-- version 2: se_cup (k+1/2) = se(k) => smoother profiles
            ! do i=its,itf
            !     if(ierr(i) /= 0) cycle
            !     do k=kts,ktf
            !         se_cup_chem(1:mtp,i,k)=se_chem(1:mtp,i,k)
            !     enddo
            ! enddo
!$acc parallel loop gang vector collapse(3)
            do k = kts, ktf
                do i = its,itf
                    do m = 1,mtp
                        if(ierr(i) == 0) se_cup_chem(m,i,k)=se_chem(m,i,k)
                    enddo
                enddo
            enddo
!$acc end parallel
        ENDIF
   
   !------------------------------------------------------------------------------------
   
    END SUBROUTINE cup_env_clev_chem

    subroutine fct1d3 (ktop,n,dt,z,tracr,massflx,trflx_in,trflx_out)

        ! --- modify a 1-D array of tracer fluxes for the purpose of maintaining
        ! --- monotonicity (including positive-definiteness) in the tracer field
        ! --- during tracer transport.
      
        ! --- the underlying transport equation is   (d tracr/dt) = - (d trflx/dz)
        ! --- where  dz = |z(k+1)-z(k)| (k=1,...,n) and  trflx = massflx * tracr
      
        ! --- note: tracr is carried in grid cells while z and fluxes are carried on
        ! --- interfaces. interface variables at index k are at grid location k-1/2.
        ! --- sign convention: mass fluxes are considered positive in +k direction.
      
        ! --- massflx and trflx_in  must be provided independently to allow the
        ! --- algorithm to generate an auxiliary low-order (diffusive) tracer flux
        ! --- as a stepping stone toward the final product trflx_out.  
      
        implicit none
        integer,intent(IN) :: n,ktop                        ! number of grid cells
        real   ,intent(IN) :: dt                            ! transport time step
        real   ,intent(IN) :: z(n+0)                        ! location of cell interfaces 
        real   ,intent(IN) :: tracr(n)                      ! the transported variable
        real   ,intent(IN) :: massflx  (n+0)                ! mass flux across interfaces
        real   ,intent(IN) :: trflx_in (n+0)                ! original tracer flux
        real   ,intent(OUT):: trflx_out(n+0)                ! modified tracr flux
        integer k,km1,kp1
        logical :: NaN, error=.false., vrbos=.false.
        real dtovdz(n),trmax(n),trmin(n),flx_lo(n+0),antifx(n+0),clipped(n+0),  &
            soln_hi(n),totlin(n),totlout(n),soln_lo(n),clipin(n),clipout(n),arg
        real,parameter :: epsil=1.e-22           ! prevent division by zero
        real,parameter :: damp=1.                ! damper of antidff flux (1=no damping)
         
        NaN(arg) = .not. (arg.ge.0. .or. arg.le.0.)        ! NaN detector
        soln_lo(:)=0.
        antifx (:)=0.
        clipout(:)=0.
         
        do k=1,ktop
            dtovdz(k)=.01*dt/abs(z(k+1)-z(k))                ! time step / grid spacing
            if (z(k).eq.z(k+1)) error=.true.
        end do
        if (vrbos .or. error) print '(a/(8es10.3))','(fct1d) dtovdz =',dtovdz(1:ktop)
         
            do k=2,ktop
                if (massflx(k).gt.0.) then
                    flx_lo(k)=massflx(k)*tracr(k-1)              ! low-order flux, upstream
                else
                    flx_lo(k)=massflx(k)*tracr(k)                ! low-order flux, upstream
                end if
                antifx(k)=trflx_in(k)-flx_lo(k)                ! antidiffusive flux
            end do
            flx_lo(  1)   =trflx_in(  1)
            flx_lo(ktop+1)=trflx_in(ktop+1)
            antifx(  1)   =0.
            antifx(ktop+1)=0.
            ! --- clip low-ord fluxes to make sure they don't violate positive-definiteness
            do k=1,ktop
                totlout(k)=max(0.,flx_lo(k+1))-min(0.,flx_lo(k  ))         ! total flux out
                clipout(k)=min(1.,tracr(k)/max(epsil,totlout(k))/ (1.0001*dtovdz(k)))
            end do
      
            do k=2,ktop
                if (massflx(k).ge.0.)  then
                    flx_lo(k)=flx_lo(k)*clipout(k-1)
                else
                    flx_lo(k)=flx_lo(k)*clipout(k)
                end if
            end do
            if (massflx(  1).lt.0.   ) flx_lo(  1)   =flx_lo(  1)*clipout(1)
            if (massflx(ktop+1).gt.0.) flx_lo(ktop+1)=flx_lo(ktop+1)*clipout(ktop)
      
            ! --- a positive-definite low-order (diffusive) solution can now be  constructed
        
            do k=1,ktop
                soln_lo  (k)=tracr(k)-(flx_lo(k+1)-flx_lo(k))*dtovdz(k)        ! low-ord solutn
                trflx_out(k)=-g*(flx_lo(k+1)-flx_lo(k))*dtovdz(k)/dt
            end do
      
            return
            soln_hi(:)=0.
            clipin (:)=0.
         
            do k=1,ktop
                km1=max(1,k-1)
                kp1=min(n,k+1)
                trmax(k)=       max(soln_lo(km1),soln_lo(k),soln_lo(kp1),        &
                                    tracr  (km1),tracr  (k),tracr  (kp1))        ! upper bound
                trmin(k)=max(0.,min(soln_lo(km1),soln_lo(k),soln_lo(kp1),        &
                                    tracr  (km1),tracr  (k),tracr  (kp1)))       ! lower bound
            end do
      
            do k=1,ktop
                totlin (k)=max(0.,antifx(k  ))-min(0.,antifx(k+1))                ! total flux in
                totlout(k)=max(0.,antifx(k+1))-min(0.,antifx(k  ))                ! total flux out
            
                clipin (k)=min(damp,(trmax(k)-soln_lo(k))/max(epsil,totlin (k))        &
                    / (1.0001*dtovdz(k)))
                clipout(k)=min(damp,(soln_lo(k)-trmin(k))/max(epsil,totlout(k))        &
                    / (1.0001*dtovdz(k)))
        
                if (NaN(clipin(k))) print *,'(fct1d) error: clipin is NaN,  k=',k
                if (NaN(clipout(k))) print *,'(fct1d) error: clipout is NaN,  k=',k
        
                if (clipin(k).lt.0.) then
                    print 100,'(fct1d) error: clipin < 0 at k =',k,                        &
                    'clipin',clipin(k),'trmax',trmax(k),'soln_lo',soln_lo(k),        &
                    'totlin',totlin(k),'dt/dz',dtovdz(k)
                    error=.true.
                end if
                if (clipout(k).lt.0.) then
                    print 100,'(fct1d) error: clipout < 0 at k =',k,                        &
                    'clipout',clipout(k),'trmin',trmin(k),'soln_lo',soln_lo(k),        &
                    'totlout',totlout(k),'dt/dz',dtovdz(k)
                    error=.true.
                end if
                100 format (a,i3/(4(a10,"=",es9.2)))
            end do
      
            do k=2,ktop
                if (antifx(k).gt.0.)  then
                    clipped(k)=antifx(k)*min(clipout(k-1),clipin(k))
                else
                    clipped(k)=antifx(k)*min(clipout(k),clipin(k-1))
                end if
                trflx_out(k)=flx_lo(k)+clipped(k)
                if (NaN(trflx_out(k)))  then
                    print *,'(fct1d) error: trflx_out is NaN,  k=',k
                    error=.true.
                end if
            end do
            do k=2,ktop
                soln_hi(k)=tracr(k)-(trflx_out(k+1)-trflx_out(k))*dtovdz(k)
                write(32,*)'3',k,soln_lo(k),soln_hi(k)
            end do
            trflx_out(  1)=trflx_in(  1)
            trflx_out(ktop+1)=trflx_in(ktop+1)
      
            if (vrbos .or. error) then
                do k=2,ktop
                    write(32,99)k,                   &
                    'tracr(k)', tracr(k),            &
                    'flx_in(k)', trflx_in(k),        &
                    'flx_in(k+1)', trflx_in(k+1),    &
                    'flx_lo(k)', flx_lo(k),          &
                    'flx_lo(k+1)', flx_lo(k+1),      &
                    'soln_lo(k)', soln_lo(k),        &
                    'trmin(k)', trmin(k),            &
                    'trmax(k)', trmax(k),            &
                    'totlin(k)', totlin(k),          &
                    'totlout(k)', totlout(k),        &
                    'clipin(k-1)', clipin(k-1),      &
                    'clipin(k)', clipin(k),          &
                    'clipout(k-1)', clipout(k-1),    &
                    'clipout(k)', clipout(k),        &
                    'antifx(k)', antifx(k),          &
                    'antifx(k+1)', antifx(k+1),      &
                    'clipped(k)', clipped(k),        &
                    'clipped(k+1)', clipped(k+1),    &
                    'flx_out(k)', trflx_out(k),      &
                    'flx_out(k+1)', trflx_out(k+1),  &
                    'dt/dz(k)', dtovdz(k),           &
                    'final', tracr(k)-(trflx_out(k+1)-trflx_out(k))*dtovdz(k)
            99    format ('(trc1d)   k =',i4/(3(a13,'=',es13.6)))
                end do
                if (error) stop '(fct1d error)'
            end if
             
        return
    end subroutine fct1d3

    SUBROUTINE bidiag (m,b,c,f)
        !-- this routine solves the problem:  bb*f(k,t+1) + cc*f(k+1,t+1) = dd
        !-- an updated "f" at time t+1 is the output
        !$acc routine seq
        implicit none
        integer, intent(in) :: m
        real, dimension(m), intent(inout) :: b,c
        real, dimension(m), intent(inout) :: f
        !--locals
        real, dimension(m) :: q
        integer :: k
        real :: p
        
        c(m)=0.
        q(1)=-c(1)/b(1)
        f(1)= f(1)/b(1)
        do k=2,m
            p	 = 1./b(k)
            q(k) = -c(k)*p
            f(k) =  f(k)*p
        enddo
        do k=m-1,1,-1
            f(k) = f(k) +q(k)*f(k+1)
        enddo
       
    END SUBROUTINE bidiag

    function auto_rk(n,step,aux,xexp,qrc1) RESULT(PW)
        !$acc routine seq
        INTEGER, INTENT(in) :: n
        real   , INTENT(in) :: step,aux,qrc1,xexp
        real                :: PW
  
        PW=step*qrc1*(1.0-exp(-aux**xexp))/float(n)
  
    end function auto_rk

    SUBROUTINE cup_env(z,qes,he,hes,t,q,p,z1,psur,ierr,itest,                     &
        itf,ktf,its,ite, kts,kte             )
            implicit none

            integer ,intent (in)                ::                         &
                itf,ktf,its,ite, kts,kte
            !
            ! ierr error value, maybe modified in this routine
            ! q           = environmental mixing ratio
            ! qes         = environmental saturation mixing ratio
            ! t           = environmental temp
            ! tv          = environmental virtual temp
            ! p           = environmental pressure
            ! z           = environmental heights
            ! he          = environmental moist static energy
            ! hes         = environmental saturation moist static energy
            ! psur        = surface pressure
            ! z1          = terrain elevation
            !
            !
            real,    dimension (its:ite,kts:kte)                              &
               ,intent (in   )                   ::                           &
                p,t,q
            real,    dimension (its:ite,kts:kte)                              &
               ,intent (out  )                   ::                           &
                he,hes,qes
            real,    dimension (its:ite,kts:kte)                              &
               ,intent (inout)                   ::                           &
               z
            real,    dimension (its:ite)                                      &
               ,intent (in   )                   ::                           &
               psur,z1
            integer, dimension (its:ite)                                      &
               ,intent (inout)                   ::                           &
               ierr
            integer                                                           &
               ,intent (in   )                   ::                           &
               itest
            !
            !  local variables in this routine
            !

            integer                              ::                           &
                i,k,iph
            !     real, dimension (1:2) :: AE,BE,HT
            real, dimension (its:ite,kts:kte) :: tv
            real :: e,tvbar
            !      real, external :: satvap
            !      real :: satvap

            he  =0.0
            hes =0.0
            qes =0.0

            !      HT(1)=xlv/cp
            !      HT(2)=2.834E6/CP
            !      BE(1)=.622*HT(1)/.286
            !      AE(1)=BE(1)/273.+ALOG(610.71)
            !      BE(2)=.622*HT(2)/.286
            !      AE(2)=BE(2)/273.+ALOG(610.71)
            !      print *, 'TCRIT = ', tcrit,its,ite
!$acc parallel
!$acc loop gang vector collapse(2) private(e)
            do k=kts,ktf
                do i=its,itf
                    if(ierr(i).eq.0)then
                        !Csgb - IPH is for phase, dependent on TCRIT (water or ice)
                        !       IPH=1
                        !       IF(T(I,K).LE.TCRIT)IPH=2
                        !       print *, 'AE(IPH),BE(IPH) = ',AE(IPH),BE(IPH),AE(IPH)-BE(IPH),T(i,k),i,k
                        !       E=EXP(AE(IPH)-BE(IPH)/T(I,K))
                        !       print *, 'P, E = ', P(I,K), E
                        !       QES(I,K)=.622*E/(100.*P(I,K)-E)
                        !
                        e=satvap(t(i,k))
                        qes(i,k)=0.622*e/max(1.e-8,(p(i,k)-e))
                        !print*,"qes=",k,qes(i,k),e,p(i,k);flush(6)
                        IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
                        IF(QES(I,K).GT.0.5   )QES(I,K)=0.5
                        IF(QES(I,K).LT.Q(I,K))QES(I,K)=Q(I,K)
                        !       IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
                        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
                    endif
                enddo
            enddo
            !
            !--- z's are calculated with changed h's and q's and t's
            !--- if itest=2
            !
            if(itest.eq.1 .or. itest.eq.0)then
!$acc loop vector
                do i=its,itf
                    if(ierr(i).eq.0)then
                        Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))- &
                        ALOG(PSUR(I)))*287.*TV(I,1)/g
                    endif
                enddo

                ! --- calculate heights
!$acc loop gang vector collapse(2) private(TVBAR)
                do K=kts+1,ktf
                    do i=its,itf
                        if(ierr(i).eq.0)then
                            TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
                            Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
                                ALOG(P(I,K-1)))*287.*TVBAR/g
                        endif
                    enddo
                enddo
            else if(itest.eq.2)then
!$acc loop gang vector collapse(2)
                do k=kts,ktf
                    do i=its,itf
                        if(ierr(i).eq.0)then
                            z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/g
                            z(i,k)=max(1.e-3,z(i,k))
                        endif
                    enddo
                enddo
            else if(itest.eq.-1)then
            endif
            !
            !--- calculate moist static energy - HE
            !    saturated moist static energy - HES
            !
!$acc loop gang vector collapse(2)
            do k=kts,ktf
                do i=its,itf
                    if(ierr(i).eq.0)then
                        if(itest.le.0)HE (I,K)=g*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
                        HES(I,K)=g*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
                        !print*,"HES=",k,HES(i,k),HE(I,K),Q(i,k),QES(I,K);flush(6)
                        IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
                    endif
                enddo
            enddo

!$acc end parallel

    END SUBROUTINE cup_env

    REAL FUNCTION satvap(temp2)
    !$acc routine seq
        implicit none
        real :: temp2, temp, toot, toto, eilog, tsot,  &
                ewlog, ewlog2, ewlog3, ewlog4
        temp = temp2-273.155
        if (temp.lt.-20.) then   !!!! ice saturation
            toot = 273.16 / temp2
            toto = 1 / toot
            eilog = -9.09718 * (toot - 1) - 3.56654 * (log(toot) / &
            log(10.)) + .876793 * (1 - toto) + (log(6.1071) / log(10.))
            satvap = 10 ** eilog
        else
            tsot = 373.16 / temp2
            ewlog = -7.90298 * (tsot - 1) + 5.02808 * &
                    (log(tsot) / log(10.))
            ewlog2 = ewlog - 1.3816e-07 * &
                    (10 ** (11.344 * (1 - (1 / tsot))) - 1)
            ewlog3 = ewlog2 + .0081328 * &
                    (10 ** (-3.49149 * (tsot - 1)) - 1)
            ewlog4 = ewlog3 + (log(1013.246) / log(10.))
            satvap = 10 ** ewlog4
        endif

    END FUNCTION
      
    SUBROUTINE get_lcl(t0,pp0,r0,tlcl,plcl,dzlcl)
        !$acc routine seq
        implicit none
        real,INTENT(IN ) :: t0,pp0,r0
        real,INTENT(OUT) :: tlcl,plcl,dzlcl
    
        real, parameter :: &
                cpg      = 102.45             &
            ,   rgas     = 287.               &
            ,   cp       = 1004.              &
            ,   p00      = 1.e5               &
            ,   g        = 9.80               &
            ,   rocp     = rgas / cp          &
            ,   p00i     = 1. / p00           &
            ,   cpor     = cp / rgas          &
            ,   cpi      = 1. / cp            &
            ,   p00k     = 26.870941          &  !  = p00 ** rocp
            ,   p00ki    = 1. / p00k
    
        integer :: nitt,ip
        real :: p0k,pi0i,ttth0,ttd,dz,pki,pppi,ti,rvs,e
        !real, external :: td,satvap
    
        !================
        !-simpler,cheaper method
        ttd=td(pp0,r0)
        tlcl=ttd-(0.001296*ttd+0.1963)*(t0-ttd)
        plcl=pp0*(tlcl/t0)**cpor
        dzlcl=127*(t0-ttd)
        if(dzlcl.le.0.)dzlcl=-999.
        !print*,"1st meth",tlcl,plcl,dzlcl;flush(6)
        RETURN
        !================
        !-2nd method
        dzlcl=-999.
        ip=0
        11 continue
    
        plcl=pp0
        tlcl=t0
        p0k=pp0**rocp
        pi0i=p0k/p00k*cp
        ttth0=t0*p00k/p0k
        ttd=td(pp0,r0)
        dz=cpg*(t0-ttd)
        if(dz.le.0.)then
            dzlcl=-999.
            return
        endif
        do 100 nitt=1,50
            pki=pi0i-g*dz/(ttth0*(1.+.61*r0))
            pppi=(pki/cp)**cpor*p00
            ti=ttth0*pki/cp
            e=100.*satvap(ti)
            rvs= ( 0.622*e )/ max(1.e-8,(pppi-e))
            !print*,'1',nitt,rvs,r0,ttd,ti,dz
            if(abs(rvs-r0).lt..00003)go to 110
            ttd=td(pppi,r0)
            dz=dz+cp/g*(ti-ttd)
            !print*,'2',nitt,rvs-r0,ttd,ti,dz
        100 continue
        print*, 'no converge in LCL:',t0,pp0,r0
        ip=ip+1
        if(ip==1)go to 11
        RETURN
    
        110 continue
        !- solution for LCL
        plcl=pppi
        tlcl=ti
        dzlcl=dz !displacement
        !print*,"2nd meth",tlcl,plcl,dz
    END SUBROUTINE get_lcl

    real function td(p,rs)
        !$acc routine seq
        implicit none
        real :: rr,rs,es,esln,p
        rr=rs+1e-8
        es=p*rr/(.622+rr)
        esln=log(es)
        td=(35.86*esln-4947.2325)/(esln-23.6837)
        return
    end function td

    SUBROUTINE cup_up_aa0(aa0,z_cup,zu,dby,GAMMA_CUP,t_cup,   &
        k22,klcl,kbcon,ktop,ierr,                      &
        itf,ktf,its,ite, kts,kte,integ_interval        )
            ! aa0 cloud work function
            ! gamma_cup = gamma on model cloud levels
            ! t_cup = temperature (Kelvin) on model cloud levels
            ! dby = buoancy term
            ! zu= normalized updraft mass flux
            ! z = heights of model levels
            ! ierr error value, maybe modified in this routine
            !
            IMPLICIT NONE
            !
            ! on input
            integer,intent (in   )              ::  itf,ktf,its,ite, kts,kte
            real,    dimension (its:ite,kts:kte) ,intent (in   )  ::  &
            z_cup,zu,gamma_cup,t_cup,dby
            integer, dimension (its:ite)         ,intent (in   )  ::  &
            k22,klcl,kbcon,ktop
            character(LEN=*),OPTIONAL, intent(in) :: integ_interval
            ! input and output
            integer, dimension (its:ite) ,intent (inout)  :: ierr
            real,    dimension (its:ite) ,intent (out  )  :: aa0
            !
            !  local variables in this routine
            integer                             ::  i,k
            real                                ::  dz,da,aa_2,aa_1
            integer, dimension (its:ite)        ::  kbeg,kend
            !
            !
            !  initialize array to zero.
            aa0(:)=0.
            !  set domain of integration
            IF(present(integ_interval)) then
                if(integ_interval == 'BL') then
                    kbeg(:) = kts
                    kend(:) = kbcon(:)-1
                elseif(integ_interval == 'CIN') then
                    kbeg(:) = k22(:) !klcl (:) ! kts
                    kend(:) = kbcon(:)-1
                else
                    STOP "unknown range in cup_up_aa0"
                endif
            ELSE
                kbeg(:) = kbcon(:)
                kend(:) = ktop (:)
            ENDIF

            do i=its,itf
                if(ierr(i) /= 0)cycle
                do k= kbeg(i),kend(i)

                    dz=z_cup(i,k+1)-z_cup(i,k)
                    aa_1=zu(i,k  )*(g/(cp*t_cup(i,k  )))*dby(i,k  )/(1.+gamma_cup(i,k  ))
                    aa_2=zu(i,k+1)*(g/(cp*t_cup(i,k+1)))*dby(i,k+1)/(1.+gamma_cup(i,k+1))
                    da=0.5*(aa_1+aa_2)*dz

                    aa0(i)=aa0(i)+da

                    !aa0(i)=aa0(i)+max(0.,da)
                enddo
            enddo

    END SUBROUTINE cup_up_aa0

    SUBROUTINE get_partition_liq_ice(ierr,tn,z1,zo_cup,po_cup, p_liq_ice,melting_layer         &
        ,itf,ktf,its,ite, kts,kte, cumulus          )
        IMPLICIT NONE
        CHARACTER *(*), INTENT (IN)                          :: cumulus
        INTEGER  ,INTENT (IN   )                             :: itf,ktf, its,ite, kts,kte
        INTEGER  ,INTENT (IN   ), DIMENSION(its:ite)         :: ierr
        REAL     ,INTENT (IN   ), DIMENSION(its:ite)         :: z1
        REAL     ,INTENT (IN   ), DIMENSION(its:ite,kts:kte) :: tn,po_cup,zo_cup
        REAL     ,INTENT (INOUT), DIMENSION(its:ite,kts:kte) :: p_liq_ice,melting_layer
        INTEGER :: i,k
        REAL    :: dp, height
        REAL, DIMENSION(its:ite) :: norm
        REAL, PARAMETER ::  T1=276.16, Z_meltlayer1=4000.,Z_meltlayer2=6000.,delT=3.
        p_liq_ice    (:,:) = 1.
        melting_layer(:,:) = 0.
        !-- get function of T for partition of total condensate into liq and ice phases.
        IF(MELT_GLAC .and. trim(cumulus) == 'deep') then
            DO k=kts,ktf
                DO i=its,itf
                    if(ierr(i) /= 0) cycle
                    if    (tn(i,k) <= T_ice) then
                        p_liq_ice(i,k) = 0.
                    elseif(  tn(i,k) > T_ice .and. tn(i,k) < T_0) then
                        p_liq_ice(i,k) =  ((tn(i,k)-T_ice)/(T_0-T_ice))**2
                    else
                        p_liq_ice(i,k) = 1.
                    endif
                    !p_liq_ice(i,k) =  min(1., (max(0.,(tn(i,k)-T_ice))/(T_0-T_ice))**2)
                ENDDO
            ENDDO
            !        go to 650
            !
            !-- define the melting layer (the layer will be between T_0+1 < TEMP < T_1
            !-- definition em terms of temperatura
            DO k=kts,ktf
                DO i=its,itf
                    if(ierr(i) /= 0) cycle
                    if    (tn(i,k) <= T_0-delT) then
                        melting_layer(i,k) = 0.

                    elseif(  tn(i,k) < T_0+delT .and. tn(i,k) > T_0-delT) then
                        melting_layer(i,k) =  ((tn(i,k)-(T_0-delt))/(2.*delT))**2

                    else
                        melting_layer(i,k) = 1.
                    endif
                        melting_layer(i,k) = melting_layer(i,k)*(1.-melting_layer(i,k))
                ENDDO
            ENDDO
            !go to 655
            !650 continue
            !        !-- definition em terms of height above local terrain
            !        DO k=kts,ktf
            !          DO i=its,itf
            !             if(ierr(i) /= 0) cycle
            !             height= zo_cup(i,k)+z1(i)
            !             if   (height > Z_meltlayer2 ) then
            !                melting_layer(i,k) = 0.
            !
            !             elseif(height > Z_meltlayer1  .and. height < Z_meltlayer2 ) then
            !
            !                melting_layer(i,k) =  ((height - Z_meltlayer1)/(Z_meltlayer2-Z_meltlayer1))**2.
            !
            !
            !             else
            !                melting_layer(i,k) = 1.
            !             endif
            !             melting_layer(i,k) = melting_layer(i,k)*(1.-melting_layer(i,k))
            !          ENDDO
            !        ENDDO
            !
            !
            !         655 continue
            !-normalize vertical integral of melting_layer to 1
            norm(:)=0.
            DO k=kts,ktf-1
                DO i=its,itf
                if(ierr(i) /= 0) cycle
                dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
                norm(i) = norm(i) + melting_layer(i,k)*dp/g
                ENDDO
            ENDDO
            DO i=its,itf
                if(ierr(i) /= 0) cycle
                melting_layer(i,:)=melting_layer(i,:)/(norm(i)+1.e-6)*(100*(po_cup(i,kts)-po_cup(i,ktf))/g)
                !print*,"i2=",i,maxval(melting_layer(i,:)),minval(melting_layer(i,:)),norm(i)
            ENDDO
            !--check
            !       norm(:)=0.
            !        DO k=kts,ktf-1
            !          DO i=its,itf
            !             dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
            !             norm(i) = norm(i) + melting_layer(i,k)*dp/g/(100*(po_cup(i,kts)-po_cup(i,ktf))/g)
            !             !print*,"n=",i,k,norm(i)
            !          ENDDO
            !        ENDDO

            !~ ELSE
            !~ p_liq_ice    (:,:) = 1.
            !~ melting_layer(:,:) = 0.
        ENDIF
    END  SUBROUTINE get_partition_liq_ice

    SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,              &
        itf,ktf,its,ite, kts,kte                     )

        IMPLICIT NONE
        !
        !  on input
        !

        ! only local dimensions are need as of now in this routine

        integer                                                           &
            ,intent (in   )                   ::                           &
            itf,ktf,                                    &
            its,ite, kts,kte
        ! array input array
        ! x output array with return values
        ! kt output array of levels
        ! ks,kend  check-range
        real,    dimension (its:ite,kts:kte)                              &
            ,intent (in   )                   ::                           &
            array
        integer, dimension (its:ite)                                      &
            ,intent (in   )                   ::                           &
            ierr,ks,kend
        integer, dimension (its:ite)                                      &
            ,intent (out  )                   ::                           &
            kt
        real,    dimension (its:ite)         ::                           &
            x
        integer                              ::                           &
            i,k,kstop

        DO 200 i=its,itf
            KT(I)=KS(I)
            if(ierr(i).eq.0)then
                X(I)=ARRAY(I,KS(I))
                KSTOP=MAX(KS(I)+1,KEND(I))
            !
                DO 100 K=KS(I)+1,KSTOP
                    IF(ARRAY(I,K).LT.X(I)) THEN
                        X(I)=ARRAY(I,K)
                        KT(I)=K
                    ENDIF
                100  CONTINUE
            endif
        200  CONTINUE

    END SUBROUTINE cup_MINIMI

    SUBROUTINE rates_up_pdf(name,ktop,ierr,p_cup,entr_rate_2d,hkbo,heo,heso_cup,z_cup, &
        kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl,klcl,hcot)

        implicit none
        integer, intent(in) :: its,ite,itf,kts,kte,ktf
        real, dimension (its:ite,kts:kte),intent (inout) :: entr_rate_2d,zuo
        real, dimension (its:ite,kts:kte),intent (in) ::p_cup, heo,heso_cup,z_cup
        real, dimension (its:ite),intent (in) :: hkbo
        integer, dimension (its:ite),intent (in) :: kstabi,k22,kbcon,kpbl,klcl
        integer, dimension (its:ite),intent (inout) :: ierr,ktop
        real, dimension (its:ite,kts:kte) :: hcot
        character *(*), intent (in) ::    name
        real :: dz,dh, dbythresh
        real :: dby(kts:kte)
        integer :: i,k,ipr,kdefi,kstart,kbegzu,kfinalzu
        integer, dimension (its:ite) :: start_level
        integer,Parameter :: FIND_KTOP_OPTION = 1 !0=original, 1=new

        dbythresh=1. !0.95  ! the range of this parameter is 0-1, higher => lower
        ! overshoting (cheque AA0 calculation)
        ! rainfall is too sensible this parameter
        ! for now, keep =1.


        dby =0.0
        hcot=0.0

        if(name == 'shallow'.or. name == 'mid')then
            dbythresh=1.0
        endif
        !         print*,"================================CUMULUS=",NAME; flush(6)
        DO i=its,itf
            kfinalzu=ktf-2
            ktop(i)=kfinalzu
            if(ierr(i).eq.0)then

                start_level(i)=kbcon(i)
                !-- hcot below kbcon
                hcot(i,kts:start_level(i))=hkbo(i)

                dz=z_cup(i,start_level(i))-z_cup(i,start_level(i)-1)
                dby(start_level(i))=(hcot(i,start_level(i))-heso_cup(i,start_level(i)))*dz

                !print*,'hco1=',start_level(I),kbcon(i),hcot(i,start_level(i))/heso_cup(i,start_level(i))

                do k=start_level(i)+1,ktf-2
                    dz=z_cup(i,k)-z_cup(i,k-1)

                    hcot(i,k)=( (1.-0.5*entr_rate_2d(i,k-1)*dz)*hcot(i,k-1)    &
                                +entr_rate_2d(i,k-1)*dz *heo (i,k-1) )/ &
                        (1.+0.5*entr_rate_2d(i,k-1)*dz)
                    dby(k)=dby(k-1)+(hcot(i,k)-heso_cup(i,k))*dz
                    !print*,'hco2=',k,hcot(i,k)/heso_cup(i,k),dby(k),entr_rate_2d(i,k-1)

                enddo
                if(FIND_KTOP_OPTION==0) then
                    do k=maxloc(dby(:),1),ktf-2
                        !~ print*,'hco30=',k,dby(k),dbythresh*maxval(dby)

                        if(dby(k).lt.dbythresh*maxval(dby))then
                            kfinalzu = k - 1
                            ktop(i)  = kfinalzu
                            !print*,'hco4=',k,kfinalzu,ktop(i),kbcon(i)+1;flush(6)
                            go to 412
                        endif
                    enddo
412                 continue
                else
                    do k=start_level(i)+1,ktf-2
                    !~ print*,'hco31=',k,dby(k),dbythresh*maxval(dby)

                        if(hcot(i,k) < heso_cup(i,k) )then
                            kfinalzu = k - 1
                            ktop(i)  = kfinalzu
                            !print*,'hco40=',k,kfinalzu,ktop(i),kbcon(i)+1;flush(6)
                            exit
                        endif
                    enddo
                endif
                if(kfinalzu.le.kbcon(i)+1) ierr(i)=41
                !~ print*,'hco5=',k,kfinalzu,ktop(i),kbcon(i)+1,ierr(i);flush(6)
                !
            ENDIF
        ENDDO
    END SUBROUTINE rates_up_pdf

    subroutine get_bouyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
        ,hc,he_cup,hes_cup,dby,z_cup)

        IMPLICIT NONE
        integer, intent (in   )   :: itf,ktf,its,ite, kts,kte
        integer, dimension (its:ite)           ,intent (in)     ::      &
                    ierr,klcl,kbcon,ktop
        real,    dimension (its:ite,kts:kte)   ,intent (in   )  ::      &
                    hc,he_cup,hes_cup,z_cup
        real,    dimension (its:ite,kts:kte)   ,intent (OUT  )  ::      &
                    dby
        INTEGER :: i,k

        do i=its,itf
            dby (i,:)=0.
            if(ierr(i) /= 0) cycle
            do k=kts,klcl(i)
                dby (i,k)=hc (i,k)-he_cup (i,k)
            enddo
            DO  k=klcl(i)+1,ktop(i)+1
                dby (i,k)=hc (i,k)-hes_cup (i,k)
            ENDDO
        ENDDO
    end subroutine get_bouyancy

    SUBROUTINE get_melting_profile(ierr,tn_cup,po_cup, p_liq_ice,melting_layer,qrco    &
        ,pwo,edto,pwdo,melting                               &
        ,itf,ktf,its,ite, kts,kte, cumulus                   )
        IMPLICIT NONE
        CHARACTER *(*), INTENT (IN)                          :: cumulus
        INTEGER  ,INTENT (IN   )                             :: itf,ktf, its,ite, kts,kte
        INTEGER  ,INTENT (IN   ), DIMENSION(its:ite)         :: ierr
        REAL     ,INTENT (IN   ), DIMENSION(its:ite)         :: edto
        REAL     ,INTENT (IN   ), DIMENSION(its:ite,kts:kte) :: tn_cup,po_cup,qrco,pwo &
                                        ,pwdo,p_liq_ice,melting_layer
        REAL     ,INTENT (INOUT), DIMENSION(its:ite,kts:kte) :: melting
        INTEGER :: i,k
        REAL    :: dp
        REAL, DIMENSION(its:ite)         :: norm,total_pwo_solid_phase
        REAL, DIMENSION(its:ite,kts:kte) :: pwo_solid_phase,pwo_eff

        IF(MELT_GLAC .and. trim(cumulus) == 'deep') then

            norm                  = 0.0
            pwo_solid_phase       = 0.0
            pwo_eff               = 0.0
            melting               = 0.0
            !-- set melting mixing ratio to zero for columns that do not have deep convection
            DO i=its,itf
                if(ierr(i) > 0) melting(i,:) = 0.
            ENDDO

            !-- now, get it for columns where deep convection is activated
            total_pwo_solid_phase(:)=0.

            DO k=kts,ktf-1
                DO i=its,itf
                    if(ierr(i) /= 0) cycle
                    dp = 100.*(po_cup(i,k)-po_cup(i,k+1))

                    !-- effective precip (after evaporation by downdraft)
                    pwo_eff(i,k) = 0.5*(pwo(i,k)+pwo(i,k+1) + edto(i)*(pwdo(i,k)+pwdo(i,k+1)))

                    !-- precipitation at solid phase(ice/snow)
                    pwo_solid_phase(i,k) = (1.-p_liq_ice(i,k))*pwo_eff(i,k)

                    !-- integrated precip at solid phase(ice/snow)
                    total_pwo_solid_phase(i) = total_pwo_solid_phase(i)+pwo_solid_phase(i,k)*dp/g
                ENDDO
            ENDDO

            DO k=kts,ktf
                DO i=its,itf
                    if(ierr(i) /= 0) cycle
                    !-- melting profile (kg/kg)
                    melting(i,k) = melting_layer(i,k)*(total_pwo_solid_phase(i)/(100*(po_cup(i,kts)-po_cup(i,ktf))/g))
                    !print*,"mel=",k,melting(i,k),pwo_solid_phase(i,k),po_cup(i,k)
                ENDDO
            ENDDO

            !-- check conservation of total solid phase precip
            !       norm(:)=0.
            !        DO k=kts,ktf-1
            !          DO i=its,itf
            !             dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
            !             norm(i) = norm(i) + melting(i,k)*dp/g
            !          ENDDO
            !        ENDDO
            !
            !       DO i=its,itf
            !         print*,"cons=",i,norm(i),total_pwo_solid_phase(i)
            !        ENDDO
            !--

        ELSE
            !-- no melting allowed in this run
            melting     (:,:) = 0.
        ENDIF
    END  SUBROUTINE get_melting_profile

    SUBROUTINE get_lateral_massflux_down(itf,ktf, its,ite, kts,kte                    &
        ,ierr,jmin,zo_cup,zdo,xzd,zd,cdd,mentrd_rate_2d      &
        ,dd_massentro,dd_massdetro ,dd_massentr, dd_massdetr &
        ,draft,mentrd_rate,dd_massentru,dd_massdetru,lambau)

        implicit none
        character *(*), intent (in) :: draft
        integer, intent(in):: itf,ktf, its,ite, kts,kte
        real,    intent(IN):: mentrd_rate
        integer, intent(in)   , dimension(its:ite)         :: ierr,jmin
        real,    intent(in)   , dimension(its:ite        ) :: lambau
        real,    intent(in)   , dimension(its:ite,kts:kte) :: zo_cup,zdo
        real,    intent(inout), dimension(its:ite,kts:kte) :: cdd,mentrd_rate_2d,xzd,zd
        real,    intent(  out), dimension(its:ite,kts:kte) :: dd_massentro, dd_massdetro&
                                        ,dd_massentr,  dd_massdetr
        real,    intent(  out), dimension(its:ite,kts:kte), OPTIONAL &
                                        :: dd_massentru, dd_massdetru
        integer ::i,ki
        real :: dzo

        do i=its,itf
            cdd(i,:)=0.
            dd_massentr (i,:)=0.
            dd_massdetr (i,:)=0.
            dd_massentro(i,:)=0.
            dd_massdetro(i,:)=0.
            if(present(dd_massentru).and.present(dd_massdetru))then
                dd_massentru(i,:)=0.
                dd_massdetru(i,:)=0.
            endif

            if(ierr(i) /= 0) cycle

            mentrd_rate_2d(i,1:jmin(i))   =mentrd_rate
            cdd           (i,1:jmin(i)-1) =mentrd_rate
            mentrd_rate_2d(i,1)=0.

            do ki=jmin(i)   ,maxloc(zdo(i,:),1),-1

                !=> from jmin to maximum value zd -> change entrainment
                dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
                dd_massdetro(i,ki)=cdd(i,ki)*dzo*zdo(i,ki+1)
                !XXX
                dd_massentro(i,ki)= zdo(i,ki)-zdo(i,ki+1)+dd_massdetro(i,ki)
                dd_massentro(i,ki)= MAX(0.,dd_massentro(i,ki))
                !-- check dd_massdetro in case of dd_massentro has been changed above
                dd_massdetro(i,ki)= dd_massentro(i,ki)-zdo(i,ki)+zdo(i,ki+1)

                !~ if(dd_massentro(i,ki).lt.0.)then
                !~ dd_massentro(i,ki)=0.
                !~ dd_massdetro(i,ki)=zdo(i,ki+1)-zdo(i,ki)
                !~ if(zdo(i,ki+1) > 0.0)&
                !~ cdd(i,ki)=dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
                !~ endif
                !~ if(zdo(i,ki+1) > 0.0)&
                !~ mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
            enddo

            do ki=maxloc(zdo(i,:),1)-1,kts,-1
                !=> from maximum value zd to surface -> change detrainment
                dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
                dd_massentro(i,ki)=mentrd_rate_2d(i,ki)*dzo*zdo(i,ki+1)
                !XXX
                dd_massdetro(i,ki) = zdo(i,ki+1)+dd_massentro(i,ki)-zdo(i,ki)
                dd_massdetro(i,ki) = MAX(0.0,dd_massdetro(i,ki))
                !-- check dd_massentro in case of dd_massdetro has been changed above
                dd_massentro(i,ki) = dd_massdetro(i,ki)+zdo(i,ki)-zdo(i,ki+1)


                !~ if(dd_massdetro(i,ki).lt.0.)then
                !~ dd_massdetro(i,ki)=0.
                !~ dd_massentro(i,ki)=zdo(i,ki)-zdo(i,ki+1)
                !~ if(zdo(i,ki+1) > 0.0)&
                !~ mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
                !~ endif
                !~ if(zdo(i,ki+1) > 0.0)&
                !~ cdd(i,ki)= dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
            enddo

            do ki=jmin(i),kts,-1
                xzd(i,ki)= zdo(i,ki)
                zd (i,ki)= zdo(i,ki)
                dd_massentr(i,ki)=dd_massentro(i,ki)
                dd_massdetr(i,ki)=dd_massdetro(i,ki)
            enddo
            if(present(dd_massentru).and.present(dd_massdetru))then
                do ki=jmin(i),kts,-1
                    dd_massentru(i,ki)=dd_massentro(i,ki)+lambau(i)*dd_massdetro(i,ki)
                    dd_massdetru(i,ki)=dd_massdetro(i,ki)+lambau(i)*dd_massdetro(i,ki)
                enddo
            endif
        enddo

    END SUBROUTINE get_lateral_massflux_down

    SUBROUTINE cup_dd_moisture(cumulus,ierrc,zd,hcd,hes_cup,qcd,qes_cup,    &
        pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,   &
        gamma_cup,pwev,pwavo,bu,qrcd, q,he,t_cup,iloop,            &
        itf,ktf,its,ite, kts,kte                              )

        IMPLICIT NONE

        character(len=*), intent(IN) :: cumulus
        integer         ,intent (in) ::  itf,ktf,its,ite, kts,kte
        ! cdd= detrainment function
        ! q = environmental q on model levels
        ! q_cup = environmental q on model cloud levels
        ! qes_cup = saturation q on model cloud levels
        ! hes_cup = saturation h on model cloud levels
        ! hcd = h in model cloud
        ! bu = buoancy term
        ! zd = normalized downdraft mass flux
        ! gamma_cup = gamma on model cloud levels
        ! mentr_rate = entrainment rate
        ! qcd = cloud q (including liquid water) after entrainment
        ! qrch = saturation q in cloud
        ! pwd = evaporate at that level
        ! pwev = total normalized integrated evaoprate (I2)
        ! entr= entrainment rate
        !
        real,    dimension (its:ite) ,intent (in  )          ::           &
        pwavo       
        real,    dimension (its:ite,kts:kte) ,intent (in   ) ::           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,                      &
        dd_massentr,dd_massdetr,gamma_cup,q,he
        integer  ,intent (in   )                             ::           &
        iloop
        integer, dimension (its:ite) ,intent (in   )         ::           &
        jmin
        integer, dimension (its:ite) ,intent (inout)         ::           &
        ierr
        real,    dimension (its:ite,kts:kte) ,intent (out  ) :: &
        qcd,qrcd,pwd
        real,    dimension (its:ite) ,intent (out  )         :: &
        pwev,bu
        character*128 :: ierrc(its:ite)
        !
        !  local variables in this routine
        !
        integer                              ::     i,k,ki
        real                                 ::     dh,dz,dqeva,denom
        !
        bu  =0.
        pwev=0.
        qcd =0.
        qrcd=0.
        pwd =0.
        !
        !
        !
        DO 100 i=its,itf
            IF(ierr(I) /= 0) cycle
            k=jmin(i)
            dz=z_cup(i,k+1)-z_cup(i,k)
            qcd(i,k)=q_cup(i,k)

            DH=HCD(I,k)-hes_cup(I,K)

            if(dh.lt.0)then
                qrcd(i,k)=(qes_cup(i,k)+(1./xlv)*(gamma_cup(i,k) /(1.+gamma_cup(i,k)))*dh)
            else
                qrcd(i,k)=qes_cup(i,k)
            endif

            pwd (i,jmin(i))= zd(i,jmin(i))*min(0.,qcd(i,k)-qrcd(i,k))
            qcd (i,k)      = qrcd(i,k)
            pwev(i)        = pwev(i)+pwd(i,jmin(i))

            bu(i)=dz*dh
            do ki=jmin(i)-1,1,-1

                dz=z_cup(i,ki+1)-z_cup(i,ki)
                denom = (zd(i,ki+1)-0.5*dd_massdetr(i,ki)+dd_massentr(i,ki))
                if( denom ==0.0 )then
                        qcd(i,ki)= qcd(i,ki+1)
                else
                    qcd(i,ki)=(qcd(i,ki+1)*zd(i,ki+1) -0.5*dd_massdetr(i,ki)*qcd(i,ki+1)+ &
                                dd_massentr(i,ki)*q(i,ki)    ) / denom
                endif
                !
                !--- to be negatively buoyant, hcd should be smaller than hes!
                !--- ideally, dh should be negative till dd hits ground, but that is not always
                !--- the case
                !
                DH=HCD(I,ki)-hes_cup(I,Ki)
                bu(i)=bu(i)+dz*dh
                QRCD(I,Ki)=qes_cup(i,ki)+(1./XLV)*(GAMMA_cup(i,ki) /(1.+GAMMA_cup(i,ki)))*DH

                dqeva=qcd(i,ki)-qrcd(i,ki)

                if(dqeva.gt.0.)then
                    dqeva=0.
                    qrcd(i,ki)=qcd(i,ki)
                endif
                pwd(i,ki)=zd(i,ki)*dqeva ! amount of liq water evaporated
                                        ! kg[water vapor]/kg[air]

                qcd(i,ki)=qrcd(i,ki)     ! water vapor in downdraft
                pwev(i)=pwev(i)+pwd(i,ki)
            enddo
            !
            if(pwev(I).eq.0.and.iloop.eq.1)then
                ierr(i)=70
                ierrc(i)="problem with buoy in cup_dd_moisture"
            endif
            if(BU(I).GE.0.and.iloop.eq.1)then
                ierr(i)=73
                ierrc(i)="problem2 with buoy in cup_dd_moisture"
            endif
            ! Ensure that precip re-evaporation does not excede total precip
            !  if(abs(pwev(i)) > pwavo(i) )then
            !   ierr(i)=77
            !   ierrc(i)="problem 3 with evap in cup_dd_moisture"
            !  endif

        100    continue!--- end loop over i

    END SUBROUTINE cup_dd_moisture

    SUBROUTINE cup_up_aa1bl(version,aa1_bl,aa1_fa,aa1,t,tn,q,qo,dtime,po_cup, &
        z_cup,zu,dby,GAMMA_CUP,t_cup,rho,                  &
        klcl,kpbl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,&
        xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux,tau_bl,tau_ecmwf,t_star,cumulus ,tn_bl,qo_bl )

        IMPLICIT NONE
        character*(*), intent(in)                         :: cumulus
        integer      , intent (in   )                     :: &
            itf,ktf,its,ite, kts,kte,version
        ! aa0 cloud work function
        ! gamma_cup = gamma on model cloud levels
        ! t_cup = temperature (Kelvin) on model cloud levels
        ! dby = buoancy term
        ! zu= normalized updraft mass flux
        ! z = heights of model levels
        ! ierr error value, maybe modified in this routine
        !
        real,    dimension (its:ite,kts:kte) ,intent (in  ) ::    &
            z_cup,zu,gamma_cup,t_cup,dby,t,tn,q,qo,po_cup,rho,tn_bl,qo_bl
        integer, dimension (its:ite)         ,intent (in  ) ::    &
            klcl,kbcon,ktop,kpbl
        real                                 ,intent(in   ) :: dtime, t_star

        real,    dimension (its:ite)         ,intent (IN  ) ::    &
            xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux       &
            ,aa1,tau_bl,tau_ecmwf
        !
        ! input and output
        integer, dimension (its:ite),intent (inout) ::  ierr
        real,    dimension (its:ite),intent (out  ) ::  aa1_bl,aa1_fa
        !  local variables in this routine
        !
        integer                              ::    i,k,iprloc
        real                                 ::    dz,da,aa_1,aa_2,tcup,da_bl,a1_bl

        !==================
        !~ integer :: ipr
        !~ REAL, parameter:: LONPR=360.-62., LATPR=-11., AG=1.
        !~ real :: xtime=0
        !~ xtime=xtime+dtime/3600.
        !~ a1_bl=0.
        !~ iprloc=-9999
        !~ if(trim(cumulus) == 'deep') then
        !~ LOCAL: do i=its,itf
        !~ !if( ( xlats(i)> LATPR-AG .and. xlats(i) <= LATPR+AG) ) then
        !~ ! print*,"XX=",i,xlons(i),xlats(i);flush(6)
        !~ !endif
        !~ if( (xlons(i)> LONPR -AG .and. xlons(i) <= LONPR +AG) .and. &
        !~ (xlats(i)> LATPR -AG .and. xlats(i) <= LATPR +AG) ) then
        !~ iprloc=i
        !~ print*,"1loc=",xlons(iprloc)-360.,xlats(iprloc),iprloc,ierr(i)
        !~ exit LOCAL
        !~ endif
        !~ ENDDO LOCAL
        !~ endif
        !~ !=================
        aa1_bl (:)=0.
        aa1_fa (:)=0.
        IF(version == 0 ) then
            do i=its,itf
                !~ !=================
                !~ if(i==iprloc .and. trim(cumulus) == 'deep')  print*,"2loc=",xlons(i)-360.,xlats(i),ierr(i)
                !~ !=================

                if(ierr(i) /= 0 ) cycle
                !~ !=================
                !~ if(i==iprloc .and. trim(cumulus) == 'deep')  print*,"3loc=",xlons(i)-360.,xlats(i),ierr(i)
                !~ !=================
                !***       do k=kts,kbcon(i)
                do k=kts,kpbl(i)
                    dz = (z_cup (i,k+1)-z_cup (i,k))*g
                    da = dz*(tn(i,k)*(1.+0.608*qo(i,k))-t(i,k)*(1.+0.608*q(i,k)))/dtime

                    !~ !=================
                    !~ if(i==iprloc .and. trim(cumulus) == 'deep')  then
                    !~ if(k==kts) write(4,*) "===",xlons(i),xlats(i),h_sfc_flux(i),le_sfc_flux(i)
                    !~ if(k==kts) write(7,*) "===",xlons(i),xlats(i),h_sfc_flux(i),le_sfc_flux(i)
                    !~ WRITE(4,*),"bl1=",k,tn(i,k)-t(i,k),(qo(i,k)-q(i,k)),da,dz/g

                    !~ a1_bl=a1_bl+dz*(tn_bl(i,k)*(1.+0.608*qo_bl(i,k))-t(i,k)*(1.+0.608*q(i,k)))/dtime
                    !~ WRITE(7,*),"bl1=",k,tn_bl(i,k)-t(i,k),(qo_bl(i,k)-q(i,k)),dz*(tn_bl(i,k)*&
                    !~ (1.+0.608*qo_bl(i,k))-t(i,k)*(1.+0.608*q(i,k)))/dtime,dz/g
                    !~ flush(4);flush(7)
                    !~ endif
                    !~ !=================

                    !--
                    !               tcup=0.5*(t_cup(i,k+1)+t_cup(i,k))
                    !             da=(da/tcup)*dtime !UNIT J/kg
                    !--
                    aa1_bl(i)=aa1_bl(i)+da ! Units : J K / (kg seg)
                enddo
                !~ !=================
                !~ !if(i==int(itf/2+0.5) .and. xland(i) <0.99.and. ztexec(i) >0.01) then
                !~ if(i==iprloc .and. trim(cumulus) == 'deep') then
                !~ write(4,*) "bl2=",xtime,kpbl(i),kbcon(i),aa1_bl(i),aa1(i)
                !~ write(4,*) "bl3=",tau_bl(i),tau_ecmwf(i),(aa1(i)-tau_bl(i)*aa1_bl(i)/t_star)
                !~ write(4,*)"================================"
                !~ write(7,*) "bl2=",xtime,kpbl(i),kbcon(i),a1_bl,aa1(i)
                !~ write(7,*) "bl3=",tau_bl(i),tau_ecmwf(i),(aa1(i)-tau_bl(i)*a1_bl/t_star)
                !~ write(7,*)"================================"
                !~ flush(4);flush(7)
                !~ endif
                !~ !=================

            enddo
        ELSEIF(version==1) then
            do i=its,itf
                if(ierr(i) /= 0 ) cycle
                do k=kts,klcl(i)
                    dz = (z_cup (i,k+1)-z_cup (i,k))
                    aa_1=(g/(cp*t_cup(i,k  )))*dby(i,k  )*zu(i,k  )
                    aa_2=(g/(cp*t_cup(i,k+1)))*dby(i,k+1)*zu(i,k+1)
                    da=0.5*(aa_1+aa_2)*dz! Units : J / kg
                    aa1_bl(i)=aa1_bl(i)+da
                    write(10,*) k,dby(i,k  ),da,aa1_bl(i),zu(i,k)
                enddo
                do k=klcl(i)+1,kbcon(i)-1
                    dz = (z_cup (i,k+1)-z_cup (i,k))
                    aa_1=(g/(cp*t_cup(i,k  )))*dby(i,k  )/(1.+gamma_cup(i,k  ))*zu(i,k  )!
                    aa_2=(g/(cp*t_cup(i,k+1)))*dby(i,k+1)/(1.+gamma_cup(i,k+1))*zu(i,k+1)!
                    da=0.5*(aa_1+aa_2)*dz! Units : J / kg
                    aa1_bl(i)=aa1_bl(i)+da
                    write(10,*)k,dby(i,k  ),da,aa1_bl(i),zu(i,k)
                enddo

            enddo
        ELSE
            stop "unknown version option in routine: cup_up_aa1bl"
        ENDIF

        return

        do i=its,itf
            if(ierr(i) /= 0)cycle
            do k= kbcon(i),ktop(i)

                dz=z_cup(i,k+1)-z_cup(i,k)
                aa_1=(g/(cp*((t_cup(i,k  )))))*dby(i,k  )/(1.+gamma_cup(i,k  ))*zu(i,k)
                aa_2=(g/(cp*((t_cup(i,k+1)))))*dby(i,k+1)/(1.+gamma_cup(i,k+1))*zu(i,k+1)
                da=0.5*(aa_1+aa_2)*dz

                aa1_fa(i)=aa1_fa(i)+da
            enddo
        enddo

    END SUBROUTINE cup_up_aa1bl

    SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
        pw,ccn,pwev,edtmax,edtmin,maxens2,edtc,psum2,psumh, &
        ccnclean,rho,aeroevap,itf,ktf,ipr,jpr,          &
        its,ite, kts,kte                     )

        IMPLICIT NONE

        integer                                                           &
        ,intent (in   )                   ::                           &
        ipr,jpr,aeroevap,itf,ktf,           &
        its,ite, kts,kte
        integer, intent (in   )              ::                           &
        maxens2
        !
        ! ierr error value, maybe modified in this routine
        !
        real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,us,vs,z,p,pw
        real,    dimension (its:ite,1:maxens2)                            &
        ,intent (out  )                   ::                           &
        edtc
        real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        edt
        real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev,ccn,psum2,psumh,edtmax,edtmin
        real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean
        integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
        integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
        !
        !  local variables in this routine
        !

        integer i,k,kk
        real    einc,pef,pefb,prezk,zkbc
        real,    dimension (its:ite)         ::                           &
        vshear,sdp,vws
        real :: prop_c,pefc,aeroadd,alpha3,beta3,rhoc
        prop_c=8. !10.386
        alpha3 = 1.9
        beta3  = -1.13
        pefc=0.

        !
        !--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
        !
        ! */ calculate an average wind shear over the depth of the cloud
        !
        edt   =0.
        vws   =0.
        sdp   =0.
        vshear=0.
        edtc  =0.

        do kk = kts,ktf-1
            do 62 i=its,itf
                IF(ierr(i).ne.0)GO TO 62
                if (kk .le. min0(ktop(i),ktf) .and. kk .ge. kbcon(i)) then
                    vws(i) = vws(i)+ &
                        (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk))) &
                    +   abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * &
                        (p(i,kk) - p(i,kk+1))
                    sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
                endif
                if (kk .eq. ktf-1)vshear(i) = 1.e3 * vws(i) / sdp(i)
            62   continue
        end do
        do i=its,itf
            IF(ierr(i).eq.0)then
                pef=(1.591-.639*VSHEAR(I)+.0953*(VSHEAR(I)**2) &
                    -.00496*(VSHEAR(I)**3))
                if(pef.gt.0.9)pef=0.9
                if(pef.lt.0.1)pef=0.1
                !
                !--- cloud base precip efficiency
                !
                zkbc=z(i,kbcon(i))*3.281e-3
                prezk=.02
                if(zkbc.gt.3.)then
                    prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc &
                    *(- 1.2569798E-2+zkbc*(4.2772E-4-zkbc*5.44E-6))))
                endif
                if(zkbc.gt.25)then
                    prezk=2.4
                endif
                pefb=1./(1.+prezk)
                !           write(11,*)pefb,prezk,zkbc
                if(pefb.gt.0.9)pefb=0.9
                if(pefb.lt.0.1)pefb=0.1
                EDT(I)=1.-.5*(pefb+pef)
                if(aeroevap.gt.1)then
                    aeroadd=(ccnclean**beta3)*((psumh(i))**(alpha3-1)) !*1.e6
                    if(i.eq.ipr)write(0,*)'edt',ccnclean,psumh(i),aeroadd
                    !              prop_c=.9/aeroadd
                    prop_c=.5*(pefb+pef)/aeroadd
                    aeroadd=(ccn(i)**beta3)*((psum2(i))**(alpha3-1)) !*1.e6
                    if(i.eq.ipr)write(0,*)'edt',ccn(i),psum2(i),aeroadd,prop_c
                    aeroadd=prop_c*aeroadd
                    pefc=aeroadd
                    if(pefc.gt.0.9)pefc=0.9
                    if(pefc.lt.0.1)pefc=0.1
                    EDT(I)=1.-pefc
                    if(aeroevap.eq.2)EDT(I)=1.-.25*(pefb+pef+2.*pefc)
                endif


                !--- edt here is 1-precipeff!
                einc=.2*edt(i)
                do k=1,maxens2
                    edtc(i,k)=edt(i)+float(k-2)*einc
                enddo
            endif
        enddo
        do i=its,itf
            IF(ierr(i).eq.0)then
                do k=1,maxens2
                    EDTC(I,K)=-EDTC(I,K)*pwav(I)/PWEV(I)
                    IF(EDTC(I,K).GT.edtmax(i))EDTC(I,K)=edtmax(i)
                    IF(EDTC(I,K).LT.edtmin(i))EDTC(I,K)=edtmin(i)
                enddo
            endif
        enddo

    END SUBROUTINE cup_dd_edt

    SUBROUTINE cup_forcing_ens_3d(itf,ktf,its,ite, kts,kte,ens4,ensdim,ichoice,maxens,maxens2,maxens3&
        ,ierr,ierr2,ierr3,k22,kbcon,ktop      &
        ,xland,aa0,aa1,xaa0,mbdt,dtime        &
        ,xf_ens,mconv,qo                      &
        ,p_cup,omeg,zd,zu,pr_ens,edt          &
        ,tau_ecmwf,aa1_bl,xf_dicycle,xk_x  )

        IMPLICIT NONE

        integer                                                           &
        ,intent (in   )                   ::                           &
        itf,ktf,its,ite, kts,kte,ens4
        integer, intent (in   )              ::                           &
        ensdim,maxens,maxens2,maxens3
        !
        ! ierr error value, maybe modified in this routine
        ! pr_ens = precipitation ensemble
        ! xf_ens = mass flux ensembles
        ! massfln = downdraft mass flux ensembles used in next timestep
        ! omeg = omega from large scale model
        ! mconv = moisture convergence from large scale model
        ! zd      = downdraft normalized mass flux
        ! zu      = updraft normalized mass flux
        ! aa0     = cloud work function without forcing effects
        ! aa1     = cloud work function with forcing effects
        ! xaa0    = cloud work function with cloud effects
        ! edt     = epsilon
        ! dir     = "storm motion"
        ! mbdt    = arbitrary numerical parameter
        ! dtime   = dt over which forcing is applied
        ! iact_gr_old = flag to tell where convection was active
        ! kbcon       = LFC of parcel from k22
        ! k22         = updraft originating level
        ! ichoice       = flag if only want one closure (usually set to zero!)
        ! name        = deep or shallow convection flag
        !
        real,    dimension (its:ite,1:ensdim)                     &
            ,intent (inout)                   ::                           &
            pr_ens
        real,    dimension (its:ite,1:ensdim)                     &
            ,intent (out  )                   ::                           &
            xf_ens
        real,    dimension (its:ite,kts:kte)                              &
            ,intent (in   )                   ::                           &
            zd,zu,p_cup,qo
        real,    dimension (its:ite,kts:kte,1:ens4)                       &
            ,intent (in   )                   ::                           &
            omeg
        real,    dimension (its:ite)                                      &
            ,intent (in   )                   ::                           &
            xaa0
        real,    dimension (its:ite)                                      &
            ,intent (in   )                   ::                           &
            aa1,edt,xland
        real,    dimension (its:ite)                                      &
            ,intent (inout)                   ::                           &
            mconv
        real,    dimension (its:ite)                                      &
            ,intent (in   )                   ::                           &
            aa0 
        real,    dimension (its:ite)                                      &
            ,intent (in   )                   ::                           &
            mbdt
        real                                                              &
            ,intent (in   )                   ::                           &
            dtime
        integer, dimension (its:ite)                                      &
            ,intent (in   )                   ::                           &
            k22,kbcon,ktop
        integer, dimension (its:ite)                                      &
            ,intent (inout)                   ::                           &
            ierr,ierr2,ierr3
        integer                                                           &
            ,intent (in   )                   ::                           &
            ichoice
        real,    intent(IN)   , dimension (its:ite) :: aa1_bl,tau_ecmwf
        real,    intent(INOUT), dimension (its:ite) :: xf_dicycle,xk_x
        !- local var
        real  :: xff_dicycle
        !
        !  local variables in this routine
        !

        real,    dimension (1:maxens3)       ::                           &
            xff_ens3
        real,    dimension (its:ite)        ::                           &
            xk
        integer                              ::                           &
            i,k,nall,n,ne,nens,nens3
        real                                 ::                           &
            a1,a_ave,xff0,xomg

        real :: betajb
        integer :: kk
        real, dimension (its:ite) :: ens_adj!,xmbmax
        !
        ens_adj(:)=1.

        !--- LARGE SCALE FORCING
        !
        DO i=its,itf
            xf_ens(i,1:16)= 0.
            IF(ierr(i) /=  0)cycle

            xff0 = (AA1(I)-AA0(I))/DTIME
            !-- default
            xff_ens3(1) = max(0.,(AA1(I)   -AA0(I))/dtime)

            xff_ens3(2) = xff_ens3(1)
            xff_ens3(3) = xff_ens3(1)
            xff_ens3(16)= xff_ens3(1)
            !
            !--- more like Brown (1979), or Frank-Cohen (199?)
            !--- omeg is in Pa/s
            xomg=0.
            kk=0
            xff_ens3(4)=0.
            do k=kbcon(i)-1,kbcon(i)+1
                !- tmp betajb=(zu(i,k)-edt(i)*zd(i,k))
                betajb=1.
                if(betajb .gt. 0.)then
                    xomg=xomg-omeg(i,k,1)/g/betajb
                    kk=kk+1
                endif
            enddo
            if(kk.gt.0)xff_ens3(4)=xomg/float(kk) ! kg[air]/m^3 * m/s
            xff_ens3(4) = max(0.0, xff_ens3(4))
            xff_ens3(5) = xff_ens3(4)
            xff_ens3(6) = xff_ens3(4)
            xff_ens3(14)= xff_ens3(4)
            !
            !--- more like Krishnamurti et al.;
            !
            !mconv(i) = 0.
            !do k=k22(i),ktop(i)
            !    mconv(i)=mconv(i)+omeg(i,k,1)*(qo(i,k+1)-qo(i,k))/g
            !enddo
            !- 2nd option (assuming that omeg(ktop)*q(ktop)<< omeg(kbcon)*q(kbcon))
            mconv(i)  = -omeg(i,kbcon(i),1)*qo(i,kbcon(i))/g ! (kg[air]/m^3)*m/s*kg[water]/kg[air]

            mconv(i)  = max(0., mconv(i))
            xff_ens3(7) = mconv(i)
            xff_ens3(8) = xff_ens3(7)
            xff_ens3(9) = xff_ens3(7)
            xff_ens3(15)= xff_ens3(7)
            !
            !---- more like  Betchold et al (2014). Note that AA1 already includes the forcings tendencies
            xff_ens3(10)= AA1(i)/tau_ecmwf(i)

            xff_ens3(11)= xff_ens3(10)
            xff_ens3(12)= xff_ens3(10)
            xff_ens3(13)= xff_ens3(10)

            !
            if(ichoice == 0)then
                if(xff0 < 0.)then
                    xff_ens3( 1)=0.
                    xff_ens3( 2)=0.
                    xff_ens3( 3)=0.
                    xff_ens3(16)=0.

                    xff_ens3(10)=0.
                    xff_ens3(11)=0.
                    xff_ens3(12)=0.
                    xff_ens3(13)=0.
                endif
            endif

            xk(i)=(XAA0(I)-(AA1(I)          ))/MBDT(I)
            if(xk(i).le.0. .and. xk(i).gt.-0.1*mbdt(i)) xk(i)=-0.1*mbdt(i)
            if(xk(i).gt.0. .and. xk(i).lt.1.e-2       ) xk(i)=1.e-2
            !
            !---  over water, enfor!e small cap for some of the closures
            !
            if(xland(i).lt.0.1)then
                if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
                    xff_ens3(1:16) = ens_adj(i)*xff_ens3(1:16)
                endif
            endif
            !
            !--- special treatment for stability closures
            !
            if(xk(i).lt.0.)then
                if(xff_ens3( 1).gt.0.)xf_ens(i, 1)=max(0.,-xff_ens3( 1)/xk(i))
                if(xff_ens3( 2).gt.0.)xf_ens(i, 2)=max(0.,-xff_ens3( 2)/xk(i))
                if(xff_ens3( 3).gt.0.)xf_ens(i, 3)=max(0.,-xff_ens3( 3)/xk(i))
                if(xff_ens3(16).gt.0.)xf_ens(i,16)=max(0.,-xff_ens3(16)/xk(i))
            else
                xff_ens3(1 )=0.
                xff_ens3(2 )=0.
                xff_ens3(3 )=0.
                xff_ens3(16)=0.
            endif

            xf_ens(i,4) =max(0.,xff_ens3(4) )
            xf_ens(i,5) =max(0.,xff_ens3(5) )
            xf_ens(i,6) =max(0.,xff_ens3(6) )
            xf_ens(i,14)=max(0.,xff_ens3(14))

            a1=max(1.e-3,pr_ens(i,7) ); xf_ens(i,7) =max(0.,xff_ens3(7)/a1)
            a1=max(1.e-3,pr_ens(i,8) ); xf_ens(i,8) =max(0.,xff_ens3(8)/a1)
            a1=max(1.e-3,pr_ens(i,9) ); xf_ens(i,9) =max(0.,xff_ens3(9)/a1)
            a1=max(1.e-3,pr_ens(i,15)); xf_ens(i,15)=max(0.,xff_ens3(15)/a1)
            !print*,"xk1=",xk(i),-xff_ens3(10), xf_ens(i,10)
            if(xk(i).lt.0.)then
                xf_ens(i,10)= max(0.,-xff_ens3(10)/xk(i))
                xf_ens(i,11)= max(0.,-xff_ens3(11)/xk(i))
                xf_ens(i,12)= max(0.,-xff_ens3(12)/xk(i))
                xf_ens(i,13)= max(0.,-xff_ens3(13)/xk(i))
            else
                xf_ens(i,10)= 0.
                xf_ens(i,11)= 0.
                xf_ens(i,12)= 0.
                xf_ens(i,13)= 0.
            endif
            !print*,"xk2=",xk(i),-xff_ens3(10), xf_ens(i,10)


            if(ichoice.ge.1)then
                xf_ens(i,1:16) =xf_ens(i,ichoice)
            endif

            !--- this for clo-X experiment
            !GO TO 33334
            !---special combination for 'ensemble closure':
            !---over the land, only applies closures 1 and 10.
            if(ichoice == 0 .and. xland(i) < 0.1)then
                xf_ens(i,1:16) =0.5*(xf_ens(i,10)+xf_ens(i,1))
            endif
            33334 continue

            !------------------------------------


        ENDDO
        !-
        !- diurnal cycle mass flux
        !-
        IF(DICYCLE==1 .or. DICYCLE==6 )THEN

            DO i=its,itf
                xf_dicycle(i) = 0.
                IF(ierr(i) /=  0)cycle

                xff_dicycle  = (AA1(i)-AA1_BL(i))/tau_ecmwf(i)
                if(xk(i).lt.0) xf_dicycle(i)= max(0.,-xff_dicycle/xk(i))
                xf_dicycle(i)= xf_ens(i,10)-xf_dicycle(i)
            ENDDO
        ELSEIF( DICYCLE==2) THEN ! for trigger function only
            DO i=its,itf
                xf_dicycle(i) = 0.
                IF(ierr(i) /=  0)cycle

                xff_ens3(1) = max(0.,-(AA1_BL(I) - AA0(I))/dtime)

                if( xff_ens3(1)   <= 0. ) then
                    xf_ens    (i,:) = 0.0
                    xf_dicycle(i)   = 0.0
                    ierr(i)=45
                    !ierrc(i)="DCclosure"
                endif
            ENDDO
        ELSEIF( DICYCLE==3) THEN

            DO i=its,itf
                xf_dicycle(i) = 0.
                IF(ierr(i) /=  0)cycle
                xff_ens3(1) = 5.*max(0.,-(AA1(I) - AA0(I))/dtime)

                if(xk(i)<0.   )xf_dicycle(i)=max(0.,-xff_ens3(1)/xk(i))
                xf_ens    (i,:) = xf_dicycle(i)
                xf_dicycle(i)   = 0.0

            ENDDO
        ELSEIF( DICYCLE==4) THEN

            DO i=its,itf
                xf_dicycle(i) = 0.
                IF(ierr(i) /=  0)cycle
                !the signal "-" is to convert from Pa/s to kg/m2/s
                if(xk_x(i) > 0.) xf_dicycle(i)= max(0., -AA1_BL(I))/xk_x(i)

                xf_ens    (i,:) = xf_dicycle(i)
                xf_dicycle(i)   = 0.0
            ENDDO
        ELSEIF( DICYCLE==5) THEN

            DO i=its,itf
                xf_dicycle(i) = 0.
                IF(ierr(i) /=  0)cycle

                xff_ens3(1) =max(0.,-(AA1_BL(I) - AA0(I))/dtime)

                if(xk(i)<0.                     )xf_dicycle(i)=max(0.,-xff_ens3(1)/xk(i))
                !              if(xk(i)<0. .and. xff_ens3(1)>0.)xf_dicycle(i)=max(0.,-xff_ens3(1)/xk(i))
                xf_ens    (i,:) = xf_dicycle(i)
                xf_dicycle(i)   = 0.0
            ENDDO
        ELSE

            xf_dicycle(:)=0.0

        ENDIF
        !--


    END SUBROUTINE cup_forcing_ens_3d

    SUBROUTINE cup_forcing_ens_3d_mid(aa0,aa1,xaa0,mbdt,dtime,ierr,&
        po_cup,ktop,k22,kbcon,kpbl,ichoice, maxens,maxens3, &
        itf,ktf,its,ite, kts,kte, &
        tau_ecmwf,aa1_bl,xf_dicycle, &
        dhdt,xff_mid,zws,hc,hco,he_cup,heo_cup)

        IMPLICIT NONE
        ! aa0     = cloud work function without forcing effects
        ! aa1     = cloud work function with forcing effects
        ! xaa0    = cloud work function with cloud effects
        ! mbdt    = arbitrary numerical parameter
        ! dtime   = dt over which forcing is applied
        ! kbcon   = LFC of parcel from k22
        ! k22     = updraft originating level
        ! ichoice = flag if only want one closure
        ! name    = deep,mid or shallow convection flag
        !
        integer  ,intent (in   )    :: itf,ktf,its,ite, kts,kte,maxens,maxens3
        integer, dimension (its:ite)          ,intent (in) ::      &
            k22,kbcon,ktop,kpbl
        real,    dimension (its:ite,kts:kte)  ,intent (in) ::      &
            po_cup,dhdt,hc,hco,he_cup,heo_cup
        real,    dimension (its:ite)          ,intent (in) ::      &
            xaa0
        real,    dimension (its:ite)          ,intent (in) ::      &
            aa1,zws,mbdt,  aa0
        real                                  ,intent (in) :: dtime
        integer                               ,intent (in) :: ichoice
        real,    dimension (its:ite)          ,intent (IN) :: aa1_bl,tau_ecmwf
        integer, dimension (its:ite)          ,intent (inout) :: ierr
        real,    dimension (its:ite)          ,intent (inout) :: xf_dicycle
        real,    dimension (its:ite,1:maxens3),intent (out)   :: xff_mid
        !
        !  local variables in this routine
        !
        real,    dimension (1:maxens) :: xk
        integer                       :: i,k
        real                          :: xff_dicycle, trash, blqe,xff_ens1,mf_ens1
        DO i=its,itf
            !-initialization
            xff_mid   (i,:)= 0.
            xf_dicycle(i)  = 0.

            if(ierr(i) /= 0)cycle

            !- Think about this:
            !xff0= (AA1(I)-AA0(I))/DTIME
            !if(xff0.lt.0.) xff_dicycle = 0.

            XK(1)=(XAA0(I)-(AA1(I)))/MBDT(i)

            if(xk(1).le.0.and.xk(1).gt.-0.1*mbdt(i)) xk(1)=-0.1*mbdt(i)
            if(xk(1).gt.0.and.xk(1).lt.1.e-2       ) xk(1)=1.e-2

            !- diurnal cycle mass flux
            if(DICYCLE==1 .or. DICYCLE==6 .or. DICYCLE==0 ) then
                !----  Betchold et al (2014)
                xff_dicycle = AA1_BL(i)/tau_ecmwf(i)
                xf_dicycle(i) = max(0.,-xff_dicycle /xk(1))
            endif
            if(DICYCLE==5 ) then
                xff_ens1 = max(0.,(AA1(I)-AA0(I))/dtime)
                mf_ens1  = 0.0
                if(XK(1).lt.0.) mf_ens1=max(0.,-xff_ens1/xk(1))
                xf_dicycle(i)=-( max(0., -(AA1_BL(I)-AA0(I))/dtime/xk(1)) - mf_ens1 )
                !xf_dicycle(i)=-( max(0., -(AA1_BL(I)-AA0(I))/dtime/xk(1)) - max(0.,-max(0.,(AA1(I)-AA0(I))/dtime)/xk(1) ))
            endif

            !- closures 3 and 4 for mid
            if(xk(1) < 0.) then
                xff_mid(i,3)=max(0., -(AA1(i)/tau_ecmwf(i))/xk(1))
                xff_mid(i,4)=xf_dicycle(i)
            endif

            !IF( ichoice == 3) xf_dicycle(i) =0. ! No dicycle for congestus
        ENDDO

        do i=its,itf
            if(ierr(i) /= 0) cycle
            !- Boundary layer quasi-equilibrium (Raymond 1995)
            if(k22(i).lt.kpbl(i)+1)then
                blqe=0.
                do k=kts,kbcon(i) !- orig formulation
                    !do k=kts,kpbl(i)
                    blqe = blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
                enddo
                !trash = max((hc (i,kbcon(i))-he_cup (i,kbcon(i))),1.e1)!- orig formulation
                trash = max((hco(i,kbcon(i))-heo_cup(i,kbcon(i))),1.e1)
                xff_mid(i,2) = max(0.,blqe/trash)
            endif

            !- W* closure (Grant,2001)
            xff_mid(i,1)=0.03*zws(i)
        enddo

        !-set xf_dicycle(i)=0 in case the closure is 4 and DICYCLE closure will be applied
        IF((DICYCLE>0) .and. ichoice == 4) xf_dicycle(:)=0.

        IF( ichoice == 5) then
            if(DICYCLE>0) then
                do i=its,itf
                    if(ierr(i) /= 0) cycle
                    !xff_mid (i,5) = 0.5*xff_mid(i,4)-xff_mid(i,2)
                    xff_mid (i,5) = 0.05*xff_mid(i,4)!-xff_mid(i,2)
                    xf_dicycle(i) = 0.
                enddo
            elseif(DICYCLE==0) then
                do i=its,itf
                    if(ierr(i) /= 0) cycle
                    !xff_mid (i,5) = 0.25*xff_mid(i,4)
                    xff_mid (i,5) = 0.05*xff_mid(i,4)
                    xf_dicycle(i) = 0.
                enddo
            endif
        ENDIF

    END SUBROUTINE cup_forcing_ens_3d_mid

    SUBROUTINE  ke_to_heating(itf,ktf,its,ite, kts,kte,ktop,ierr &
        ,po_cup,us,vs,dellu,dellv,dellat)

        implicit none
        integer                              ,intent (in   ) :: itf,ktf,its,ite, kts,kte
        integer, dimension (its:ite)         ,intent (in   ) :: ierr,ktop
        real   , dimension (its:ite,kts:kte) ,intent (in   ) :: po_cup,us,vs,dellu,dellv
        real   , dimension (its:ite,kts:kte) ,intent (inout) :: dellat

        real :: dts,fp,dp,fpi
        integer ::i,k

        ! since kinetic energy is being dissipated, add heating accordingly (from ECMWF)
        !
        do i=its,itf
            if(ierr(i) /= 0) cycle
            dts=0.
            fpi=0.
            do k=kts,ktop(i)
                dp=(po_cup(i,k)-po_cup(i,k+1))*100.
                !total KE dissiptaion estimate
                dts = dts - (dellu(i,k)*us(i,k)+dellv(i,k)*vs(i,k))*dp/g
                !
                ! fpi needed for calcualtion of conversion to pot. energyintegrated
                fpi = fpi + sqrt(dellu(i,k)*dellu(i,k) + dellv(i,k)*dellv(i,k))*dp
            enddo
            if(fpi.gt.0.)then
                do k=kts,ktop(i)
                    fp= sqrt((dellu(i,k)*dellu(i,k)+dellv(i,k)*dellv(i,k)))/fpi
                    dellat(i,k) = dellat(i,k) + fp*dts*g/cp
                enddo
            endif
        enddo

    end SUBROUTINE  ke_to_heating

    SUBROUTINE cup_output_ens_3d(name,xff_mid,xf_ens,ierr,dellat,dellaq,dellaqc,  &
        outtem,outq,outqc,     &
        zu,pre,pw,xmb,ktop,                 &
        nx,nx2,ierr2,ierr3,pr_ens,             &
        maxens3,ensdim,sig,xland1,        &
        ichoice,ipr,jpr,itf,ktf,its,ite, kts,kte,   &
        xf_dicycle,outu,outv,dellu,dellv,dtime,p_cup,kbcon )

        IMPLICIT NONE
        !
        !  on input
        !

        ! only local dimensions are need as of now in this routine

        integer                                                           &
            ,intent (in   )                   ::                           &
            ichoice,ipr,jpr,itf,ktf,     &
            its,ite, kts,kte
        integer, intent (in   )              ::                           &
            ensdim,nx,nx2,maxens3
        ! xf_ens = ensemble mass fluxes
        ! pr_ens = precipitation ensembles
        ! dellat = change of temperature per unit mass flux of cloud ensemble
        ! dellaq = change of q per unit mass flux of cloud ensemble
        ! dellaqc = change of qc per unit mass flux of cloud ensemble
        ! outtem = output temp tendency (per s)
        ! outq   = output q tendency (per s)
        ! outqc  = output qc tendency (per s)
        ! pre    = output precip
        ! xmb    = total base mass flux
        ! xfac1  = correction factor
        ! pw = pw -epsilon*pd (ensemble dependent)
        ! ierr error value, maybe modified in this routine
        !
        real,    dimension (its:ite,1:ensdim)                     &
            ,intent (inout)                   ::                           &
            xf_ens,pr_ens
        real,    dimension (its:ite,kts:kte)                              &
            ,intent (out  )                   ::                           &
            outtem,outq,outqc,outu,outv
        real,    dimension (its:ite,kts:kte)                              &
            ,intent (in  )                   ::                           &
            zu,p_cup
        real,   dimension (its:ite)                                      &
            ,intent (in  )                   ::                           &
            sig
        real,   dimension (its:ite,maxens3)                                      &
            ,intent (in  )                   ::                           &
            xff_mid
        real,    dimension (its:ite)                                      &
            ,intent (out  )                   ::                           &
            pre,xmb
        real,    dimension (its:ite)                                      &
            ,intent (inout  )                   ::                           &
            xland1
        real,    dimension (its:ite,kts:kte)                     &
            ,intent (in   )                   ::                           &
            dellat,dellaqc,dellaq,pw,dellu,dellv
        integer, dimension (its:ite)                                      &
            ,intent (in   )                   ::                           &
            ktop,kbcon
        integer, dimension (its:ite)                                      &
            ,intent (inout)                   ::                           &
            ierr,ierr2,ierr3
        real,    intent(IN), dimension (its:ite) :: xf_dicycle
        real,    intent(IN) :: dtime
        !
        !  local variables in this routine
        !

        integer                              ::                           &
            i,k,n,ncount,zmax
        real                                 ::                           &
            outtes,ddtes,dtt,dtq,dtqc,dtpw,prerate,fixouts
        real                                 ::                           &
            dtts,dtqs
        real,    dimension (its:ite)::                           &
            xmb_ave,xmbmax
        character *(*), intent (in)        ::         name
        !
        DO k=kts,ktf
            do i=its,itf
                outtem (i,k)=0.
                outq   (i,k)=0.
                outqc  (i,k)=0.
                outu   (i,k)=0.
                outv   (i,k)=0.
            enddo
        enddo
        do i=its,itf
            pre(i)  =0.
            xmb(i)  =0.
        enddo

        do i=its,itf
            IF(ierr(i).eq.0)then
                do n=1,maxens3
                    if(pr_ens(i,n).le.0.)then
                        xf_ens(i,n)=0.
                    endif
                enddo
            endif
        enddo
        !
        !--- calculate ensemble average mass fluxes
        !
        IF(NAME == 'deep') THEN
            do i=its,itf
                if(ierr(i).eq.0)then
                    k=0
                    xmb_ave(i)=0.
                    do n=1,maxens3
                        k=k+1
                        xmb_ave(i)=xmb_ave(i)+xf_ens(i,n)
                    enddo
                    !- 'ensemble' average mass flux
                    xmb_ave(i)=xmb_ave(i)/float(k)
                endif
            ENDDO

            !- mid (congestus type) convection
        ELSEIF(NAME=='mid') then
            if(ichoice .le. 5) then
                do i=its,itf
                    if(ierr(i) /= 0) cycle
                    if(ichoice == 0) then
                        xmb_ave(i)=0.3333*(xff_mid(i,1)+xff_mid(i,2)+xff_mid(i,3))
                    else
                        xmb_ave(i)= xff_mid(i,ichoice)
                    endif
                enddo
            else
                stop 'For mid ichoice must be 0,..,5'
            endif
        ENDIF
        !- set the updradt mass flux and do not allow negative values and apply the diurnal cycle closure
        do i=its,itf
            if(ierr(i) /= 0) cycle
            !- mass flux of updradt at cloud base
            xmb(i) = xmb_ave(i)
            !- diurnal cycle closure
            xmb(i) = xmb(i) - xf_dicycle(i)
            if(xmb(i) .le. 0.)then
                ierr(i)=13; xmb (i)=0.
            endif
        enddo
        !-apply the scale-dependence Arakawa's approach
        do i=its,itf
            if(ierr(i) /= 0) cycle
            !- scale dependence
            xmb(i)=sig(i)*xmb(i)

            if(xmb(i) == 0. ) ierr(i)=14
            if(xmb(i) > 100.) ierr(i)=15
        ENDDO

        !--- sanity check for mass flux
        !
        DO i=its,itf
            if(ierr(i) /= 0) cycle
            xmbmax(i)=100.*(p_cup(i,kbcon(i))-p_cup(i,kbcon(i)+1))/(g*dtime)
            xmb(i) = min(xmb(i),xmbmax(i))
        enddo

        !--- check outtem and and outq for high values
        !--- criteria: if abs (dT/dt or dQ/dt) > 100 K/day => fix xmb
        DO i=its,itf
            IF(ierr(i) /= 0) CYCLE
            fixouts=xmb(i) *86400.*max(maxval(abs(dellat(i,kts:ktop(i)))),&
                                (xlv/cp)*maxval(abs(dellaq(i,kts:ktop(i)))) )

            if(fixouts > 100.) then ! K/day
                fixouts=100./(fixouts)
                xmb   (i)  = xmb   (i)  *fixouts
                xf_ens(i,:)= xf_ens(i,:)*fixouts
            endif
        ENDDO
        !
        !-- now do feedback
        !
        DO i=its,itf
            IF(ierr(i) /= 0) CYCLE
            DO k=kts,ktop(i)
                pre    (i)  = pre(i)+pw(i,k)*xmb(i)

                outtem (i,k)= dellat   (i,k)*xmb(i)
                outq   (i,k)= dellaq   (i,k)*xmb(i)
                outqc  (i,k)= dellaqc  (i,k)*xmb(i)
                outu   (i,k)= dellu    (i,k)*xmb(i)
                outv   (i,k)= dellv    (i,k)*xmb(i)

            ENDDO
            xf_ens (i,:)= sig(i)*xf_ens(i,:)
        ENDDO

    END SUBROUTINE cup_output_ens_3d

    SUBROUTINE get_precip_fluxes(cumulus,klcl,kbcon,ktop,k22,ierr,xland,pre,xmb  &
        ,pwo,pwavo,edto,pwevo,pwdo,t_cup,tempco                   &
        ,prec_flx,evap_flx                                        &
        ,itf,ktf,its,ite, kts,kte)

        implicit none
        character *(*)            , intent (in) :: cumulus        
        integer                    ,intent (in) :: itf,ktf,its,ite,kts,kte
        integer, dimension(its:ite),intent (in) :: kbcon,ktop,k22,klcl,ierr
        real,    dimension(its:ite),intent (in) :: xland,pwavo,pwevo,edto,pre,xmb
        real,    dimension(its:ite,kts:kte),intent (in)  :: pwo,pwdo,t_cup,tempco 
        real,    dimension(its:ite,kts:kte),intent (out) :: prec_flx,evap_flx !-- units kg[water]/m2/s

        !-- locals
        integer :: i,k
        prec_flx = 0.0
        evap_flx = 0.0
        if(trim(cumulus) == 'shallow') return

        DO i=its,itf
            if (ierr(i)  /= 0) cycle

            do k=ktop(i),kts,-1

                !--- precipitation flux (at 'cup' levels), units: kg[water]/m2/s
                prec_flx(i,k) = prec_flx(i,k+1) + xmb(i)*(pwo(i,k) + edto(i)*pwdo(i,k))
                prec_flx(i,k) = max(0.,prec_flx(i,k))

                !--- evaporation flux (at 'cup' levels), units: kg[water]/m2/s
                evap_flx(i,k) = evap_flx(i,k+1) - xmb(i)*edto(i)*pwdo(i,k) 
                evap_flx(i,k) = max(0.,evap_flx(i,k))


                !
                !--for future use (rain and snow precipitation fluxes)
                !if	 (tempco(i,k) <= T_ice) then
                !   p_liq_ice(k) = 0.      ! only solid 
                !elseif(  tempco(i,k) > T_ice .and. tempco(i,k) < T_0) then
                !   p_liq_ice(k) =  ((tempco(i,k)-T_ice)/(T_0-T_ice))**2
                !else
                !   p_liq_ice(k) = 1.      ! only liquid 
                !endif
                !prec_flx_rain(k) = prec_flx(i,k)*(1.-p_liq_ice(k))
                !prec_flx_snow(k) = prec_flx(i,k)*    p_liq_ice(k)
            enddo

            !if(prec_flx   (i,kts) .ne. pre(i)) then 
            !print*,"error=",100.*(prec_flx   (i,kts) - pre(i))/(1.e-16+pre(i)),pre(i),prec_flx   (i,kts)
            !STOP 'problem with water balance' 
            !endif
        ENDDO  
    END   SUBROUTINE get_precip_fluxes

    SUBROUTINE tridiag (m,a,b,c,f)
        !-- this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
        !-- an updated "f" at time t+1 is the output
        implicit none
        integer, intent(in) :: m
        real, dimension(m), intent(inout) :: a,b,c
        real, dimension(m), intent(inout) :: f
        !--locals
        real, dimension(m) :: q
        integer :: k
        real :: p
        
        c(m)=0.
        q(1)=-c(1)/b(1)
        f(1)= f(1)/b(1)
        do k=2,m
            p	 = 1./( b(k)+a(k)*q(k-1) )
            q(k) = -c(k)*p
            f(k) = p*(f(k) - a(k)*f(k-1))
        enddo
        do k=m-1,1,-1
            f(k) = f(k) +q(k)*f(k+1)
        enddo
       
    END SUBROUTINE tridiag
  
    SUBROUTINE CUP_gf_sh(itf,ktf ,its,ite, kts,kte, mtp &
        ,cumulus                            &
        ,ichoice                            &
        ,entr_rate                          &
        ,use_excess                         &
        !input data
        ,dtime                              &
        ,dx                                 &
        ,kpbl                               &
        ,h_sfc_flux                         &
        ,le_sfc_flux                        &
        ,tsur                               &
        ,psur                               &
        ,z1                                 &
        ,xland                              &
        ,ztexec                             &
        ,zqexec                             &
        ,ccn                                &
        ,rho                                &
        ,dhdt                               &
        ,zws                                &
        ,zo                                 &
        ,dm2d                               &
        ,t                                  &
        ,q                                  &
        ,tn                                 &
        ,qo                                 &
        ,po                                 &
        ,us                                 &
        ,vs                                 &
        ,se_chem                            &
        !-----output data
        ,outt                               &
        ,outq                               &
        ,outqc                              &
        ,outu                               &
        ,outv                               &
        ,out_chem                           &
       !- for convective transport
        ,ierr                               &
        ,jmin                               &
        ,klcl                               &
        ,k22                                &
        ,kbcon                              &
        ,ktop                               &
        ,kstabi                             &
        ,kstabm                             &
        ,xmb                                &
        ,edto                               &
        ,pwavo                              &
        ,po_cup                             &
        ,up_massentro                       &
        ,up_massdetro                       &
        ,dd_massentro                       &
        ,dd_massdetro                       &
        ,zuo                                &
        ,zdo                                &
        ,pwo                                &
        ,pwdo                               &
        ,qrco                               &
        ,tup                                &
        ,clfrac                             &
        )
        implicit none
        ! Note : nvfortran compiler with OpenACC does not know how to handle assumed-sized character arrays
        ! character*(*), intent(in) :: cumulus
        character*(10), intent(in) :: cumulus

        integer,intent (in   )    :: itf,ktf,its,ite, kts,kte,     &
                                    ichoice,use_excess
        !
        ! outtem = output temp tendency (per s)
        ! outq   = output q tendency (per s)
        ! outqc  = output qc tendency (per s)
        ! pre    = output precip
        real, dimension (its:ite,kts:kte),intent (inout) ::   outt,outq,outqc,outu,outv

        integer, dimension (its:ite) ,intent (in  )         ::   kpbl
        !
        ! basic environmental input includes moisture convergence (mconv)
        ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
        ! convection for this call only and at that particular gridpoint
        !
        real, dimension (its:ite,kts:kte) ,intent (in   ) ::      &
                        rho,t,po,us,vs,tn,dhdt,dm2d
        real, dimension (its:ite,kts:kte) ,intent (inout) ::      &
                        Q,QO

        real, dimension (its:ite) ,intent (in   )            ::      &
                        zws,ztexec,zqexec,ccn,z1,psur,xland        &
                        ,h_sfc_flux,le_sfc_flux,tsur,dx

        real, intent (in   )                             ::      &
                        dtime,entr_rate
        !
        !
        !***************** the following are your basic environmental
        !                     variables. They carry a "_cup" if they are
        !                     on model cloud levels (staggered). They carry
        !                     an "o"-ending (z becomes zo), if they are the forced
        !                     variables. They are preceded by x (z becomes xz)
        !                     to indicate modification by some typ of cloud
        !
        ! z           = heights of model levels
        ! q           = environmental mixing ratio
        ! qes         = environmental saturation mixing ratio
        ! t           = environmental temp
        ! p           = environmental pressure
        ! he          = environmental moist static energy
        ! hes         = environmental saturation moist static energy
        ! z_cup       = heights of model cloud levels
        ! q_cup       = environmental q on model cloud levels
        ! qes_cup     = saturation q on model cloud levels
        ! t_cup       = temperature (Kelvin) on model cloud levels
        ! p_cup       = environmental pressure
        ! he_cup = moist static energy on model cloud levels
        ! hes_cup = saturation moist static energy on model cloud levels
        ! gamma_cup = gamma on model cloud levels
        !
        !
        ! hcd = moist static energy in downdraft
        ! zd normalized downdraft mass flux
        ! dby = buoancy term
        ! entr = entrainment rate
        ! zd   = downdraft normalized mass flux
        ! entr= entrainment rate
        ! hcd = h in model cloud
        ! bu = buoancy term
        ! zd = normalized downdraft mass flux
        ! gamma_cup = gamma on model cloud levels
        ! qcd = cloud q (including liquid water) after entrainment
        ! qrch = saturation q in cloud
        ! pwd = evaporate at that level
        ! pwev = total normalized integrated evaoprate (I2)
        ! entr= entrainment rate
        ! z1 = terrain elevation
        ! entr = downdraft entrainment rate
        ! jmin = downdraft originating level
        ! kdet = level above ground where downdraft start detraining
        ! psur        = surface pressure
        ! z1          = terrain elevation
        ! pr_ens = precipitation ensemble
        ! xf_ens = mass flux ensembles
        ! massfln = downdraft mass flux ensembles used in next timestep
        ! omeg = omega from large scale model
        ! mconv = moisture convergence from large scale model
        ! zd      = downdraft normalized mass flux
        ! zu      = updraft normalized mass flux
        ! dir     = "storm motion"
        ! mbdt    = arbitrary numerical parameter
        ! dtime   = dt over which forcing is applied
        ! iact_gr_old = flag to tell where convection was active
        ! kbcon       = LFC of parcel from k22
        ! k22         = updraft originating level
        ! ichoice       = flag if only want one closure (usually set to zero!)
        ! dby = buoancy term
        ! ktop = cloud top (output)
        ! xmb    = total base mass flux
        ! hc = cloud moist static energy
        ! hkb = moist static energy at originating level

        real,    dimension (its:ite,kts:kte) ::                              &
            entr_rate_2d, tempco, he,hes,qes,z, heo,heso,qeso,zo,            &
            qes_cup,  q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup,t_cup, &
            qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,       gammao_cup,tn_cup,&
            dby,pw,hc,zu,zd,clw_all,dbyo,qco,hco,                            &

            xhe,xhes,xqes,xz,xt,xq,  xdby,xhc,xzu,                           &
            xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup, xt_cup, cupclw,         &

            ! cd  = detrainment function for updraft
            ! cdd = detrainment function for downdraft
            ! dellat = change of temperature per unit mass flux of cloud ensemble
            ! dellaq = change of q per unit mass flux of cloud ensemble
            ! dellaqc = change of qc per unit mass flux of cloud ensemble

            cd,dellah,dellaq,dellat,dellaqc

        ! aa0 cloud work function for downdraft
        ! edt = epsilon
        ! aa0     = cloud work function without forcing effects
        ! aa1     = cloud work function with forcing effects
        ! xaa0    = cloud work function with cloud effects
        ! edt     = epsilon
        integer, dimension (its:ite) :: kbmax

        real,    dimension (its:ite) :: xland1 , cap_max_increment,cap_max &
                                    ,hkb,hkbo,xhkb, xmbmax     &
                                    ,psum,psumh                &
                                    ,aa1,aa0,xaa0,cape

        integer  ::    ipr=0,jpr=0                       &
        ,nall,iedt,nens,nens3,ki,i,k,kk,iresult,start_k22 &
        ,icount,jprnt,k1,k2,nvar
        real     :: day,dz,dzo,mbdt,radius, depth_min,zkbmax,zktop,      &
        dh,cap_maxs,trash,frh,fsum, detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
        detup,subdown,entdoj,entupk,detupk,totmas,blqe,xkshal,  &
        C_up, tcold,thot,efic,fin,trash2,x_add,DELLAH_aver,    &
        zlll,tlll,plll,rlll,tlcl,plcl,dzlcl,fixouts,denomU

        real xff_shal(9)

        character*128 :: ierrc(its:ite)

        real,    dimension (its:ite,kts:kte) :: up_massentr ,up_massdetr ,dd_massentr ,dd_massdetr  &
                                            ,dtempdz,hcot,vvel2d                       &
                                            ,u_cup,v_cup,uc,vc,dellu,dellv

        real,    dimension (kts:kte) :: dummy1,dummy2

        integer, dimension (its:ite,kts:kte) ::  k_inv_layers
        real,    dimension (its:ite,kts:kte) ::  c1d
        real,    dimension (its:ite        ) ::  lambau, vvel1d
        real,    dimension (its:ite,kts:kte) ::  up_massentru,up_massdetru,&
                                                dd_massentru,dd_massdetru

        !- atmos composition treatment
        integer,intent (in   )    :: mtp
        real, dimension (mtp,its:ite,kts:kte),intent (in)    ::   se_chem
        real, dimension (mtp,its:ite,kts:kte),intent (inout) ::   out_chem

        !-locals
        real, dimension (kts:kte)     :: aa,bb,cc,ddu,ddv,ddh,ddq,fp,fm
        real :: alp0,beta1
        integer :: ispc,zmax
        real, dimension (mtp,its:ite,kts:kte) :: se_cup_chem,sc_up_chem
        real, dimension (mtp,its:ite)         :: sc_up_chem_b
        real, dimension (mtp,kts:kte)         :: trcflx_in,sub_tend,ddtr
        real, dimension (mtp)     ::   min_tend_chem,delS_up,env_sub,outchem1
        real, dimension (kts:kte) ::   distance
        real, dimension (its:ite,kts:kte)     :: massflx,zenv
        real :: XZZ,XZD,XZE,denom,dtime_max,massi,massf,env_mf
        integer, dimension (its:ite)          :: start_level

        integer, dimension (its:ite)      , intent (inout)  :: &
                ierr              &
                ,jmin              &
                ,klcl              &
                ,k22               &
                ,kbcon             &
                ,ktop              &
                ,kstabi            &
                ,kstabm

        real,  dimension (its:ite)        , intent (inout)  :: &
                xmb               &
                ,edto              &
                ,pwavo
        real,  dimension (its:ite,kts:kte), intent (inout)  :: &
                po_cup            &
                ,up_massentro      &
                ,up_massdetro      &
                ,dd_massentro      &
                ,dd_massdetro      &
                ,zuo               &
                ,zdo               &
                ,pwo               &
                ,pwdo              &
                ,qrco              &
                ,tup               &
                ,clfrac

        character(len=1):: cty
        logical, parameter :: USE_LCL       =.FALSE.
        ! logical, parameter :: USE_INV_LAYERS=.TRUE. 
        logical, parameter :: USE_INV_LAYERS=.FALSE. ! P9_P6 

        character(12) :: cumulus_1, cumulus_2

!$acc data copy(outt,outq,outqc,outu,outv, Q,QO,out_chem,ierr,jmin,klcl,k22,kbcon, &
!$acc           ktop,kstabi,kstabm,xmb,edto,pwavo,po_cup,up_massentro,up_massdetro, &
!$acc           dd_massentro,dd_massdetro,zuo,zdo,pwo,pwdo,qrco,tup,clfrac, zo) &
!$acc      copyin(kpbl, rho,t,po,us,vs,tn,dhdt,dm2d, zws,ztexec,zqexec,ccn,z1,&
!$acc             psur,xland,h_sfc_flux,le_sfc_flux,tsur,dx,se_chem) & ! Note : Cumulus does not copy in
!$acc      create(entr_rate_2d, tempco, he,hes,qes,z, heo,heso,qeso,qes_cup, &
!$acc            q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup, qeso_cup, &
!$acc            qo_cup,heo_cup,heso_cup,zo_cup,gammao_cup,tn_cup,dby,pw,hc,zu,zd, &
!$acc            clw_all,dbyo,qco,hco,xhe,xhes,xqes,xz,xt,xq,xdby,xhc,xzu,xqes_cup, &
!$acc            xq_cup,xhe_cup,xhes_cup,xz_cup, xt_cup, cupclw,cd,dellah,dellaq,&
!$acc            dellat,dellaqc,kbmax,xland1,cap_max_increment,cap_max,hkb,hkbo,xhkb, &
!$acc            xmbmax,psum,psumh,aa1,aa0,xaa0,cape,up_massentr,up_massdetr,dd_massentr,&
!$acc            dd_massdetr,dtempdz,hcot,vvel2d,u_cup,v_cup,uc,vc,dellu,dellv, &
!$acc            k_inv_layers,c1d,lambau,vvel1d,up_massentru,up_massdetru,dd_massentru,&
!$acc            dd_massdetru,aa,bb,cc,ddu,ddv,ddh,ddq,fp,fm,se_cup_chem,sc_up_chem, &
!$acc            sc_up_chem_b,trcflx_in,sub_tend,ddtr,min_tend_chem,delS_up,env_sub, &
!$acc            outchem1,distance,massflx,zenv,start_level,xff_shal, ierrc, dummy1,dummy2)

!$acc parallel loop vector
        do i=its,itf
            ierr(i)=0
            xland1(i)=xland(i) ! 1.
            if(xland(i).gt.1.5) xland1(i)=0.
            ierrc(i)=" "
            ! hcot(i,:)=0.0
            lambau(i)=2.0
        enddo
!$acc end parallel

!$acc parallel loop gang vector collapse(2)
        do k=kts,ktf
            do i=its,itf
                hcot(i,k) = 0.
                z  (i,k) = zo(i,k)
                xz (i,k) = zo(i,k)
                c1d(i,k) = 0.
            enddo
        enddo
!$acc end parallel

        !
        !
        !--- maximum depth (mb) of capping
        !--- inversion (larger cap = no convection)
        !
        !cap_maxs=150. 
        cap_maxs= 50. !P9_P6

!$acc parallel loop vector
        do i=its,itf
            cap_max(i)=cap_maxs
            cap_max_increment(i)=25.
        enddo
!$acc end parallel
        !
        !--- calculate moist static energy, heights, qes
        !
        call cup_env(z,qes,he,hes,t,q,po,z1, &
            psur,ierr,-1,   &
            itf,ktf,its,ite, kts,kte)
        call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, &
            psur,ierr,-1,   &
            itf,ktf,its,ite, kts,kte)

        !
        !--- environmental values on cloud levels
        !
        call cup_env_clev(t,qes,q,he,hes,z,po,qes_cup,q_cup,he_cup, &
            us,vs,u_cup,v_cup,                                     &
            hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,tsur,         &
            ierr,z1,itf,ktf,its,ite, kts,kte, host_model)
        call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup,heo_cup,    &
            us,vs,u_cup,v_cup,                                                 &
            heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur, tsur,               &
            ierr,z1,itf,ktf,its,ite, kts,kte, host_model)

        !
        !--- max height(m) above ground where updraft air can originate
        !
        zkbmax=3000.
!$acc kernels present(kbmax, zo_cup, z1)
        do i=its,itf
            kbmax(i)=1
            if(ierr(i).eq.0) then
                do k=kts,ktf
                    if(zo_cup(i,k).gt.zkbmax+z1(i)) then
                        kbmax(i)=k
                        go to 25
                    endif
                enddo
                25     continue

                kbmax(i)=min(kbmax(i),ktf-4)
            endif
        enddo
!$acc end kernels

        !
        !------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
        !
        start_k22=1
!$acc kernels
        k22(:)=0
!$acc end kernels

!$acc parallel loop vector
        DO 36 i=its,itf
            IF(ierr(I).eq.0) THEN

                k22(i)=maxloc(HEO_CUP(i,start_k22:kbmax(i)+2),1)+start_k22-1
                k22(i)=min(2,k22(i))

                if(K22(I).GT.KBMAX(i)) then
                    ierr(i)=2
                    ierrc(i)="could not find k22"
                endif
            endif
        36   CONTINUE
!$acc end parallel
        !
        !
        !------- DETERMINE LCL for the air parcels around K22
        !

!$acc parallel loop seq private(rlll, tlll, plll, zlll, k, x_add, tlcl, plcl, dzlcl, &
!$acc                            aa, bb, cc, ddu, ddv)
        DO i=its,itf
            klcl(i) = k22(i) ! default value

            if(ierr(i) ==0)then
                !$acc loop vector
                do k = kts,kte
                    aa(k) = po(i,k)
                    bb(k) = q_cup(i, k)
                    cc(k) = t_cup(i,k)
                    ddu(k) = p_cup(i,k)
                    ddv(k) = z_cup(i,k)
                enddo
                !- tlll, rlll,plll - temp, water vapor and pressure of the source air parcel
                x_add = float(use_excess)*zqexec(i)
                call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),aa,bb,rlll,k22(i),x_add)
                x_add = float(use_excess)*ztexec(i)
                call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),aa,cc,tlll,k22(i),x_add)
                call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),aa,ddu,plll,k22(i))

                call get_lcl(tlll,100.*plll,rlll,tlcl,plcl,dzlcl)

                !-get LCL layer index
                if(dzlcl >= 0.) then ! LCL found (if dzlcl<0 => not found)

                    call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),aa,ddv,zlll,k22(i))

                    loop0: do k=kts,ktf
                        if(z_cup(i,k).gt.zlll+dzlcl) then
                            klcl(i)=k
                            exit loop0
                        endif
                    enddo loop0
                    klcl(i)=min(klcl(i),ktf-4)
                endif
                !write(9,110)'SHlcl',tlcl,plcl,dzlcl,klcl(i)
                !110 format(1x,A5,3F10.2,i4)
            endif
        Enddo
!$acc end parallel loop

        !-- check if LCL is below PBL height for shallow convection
        IF(USE_LCL) then
            do i=its,itf
                if(ierr(i)/= 0) cycle
                if(klcl(i) > max(1,kpbl(i)+1)) then
                    ierr(i)=21
                    ierrc(i)='for shallow conv:LCL height lt PBL height'
                endif
            ENDDO
        ENDIF
        !
        !--- DETERMINE THE BOUNDARY CONDITION FOR in-cloud MOIST STATIC ENERGY
        !

!$acc kernels
        do i=its,itf
            if(ierr(I) /= 0)Cycle
            x_add = float(use_excess)*(xlv*zqexec(i)+cp*ztexec(i))
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),he_cup (i,kts:kte),hkb (i),k22(i),x_add)
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup(i,kts:kte),hkbo(i),k22(i),x_add)
        enddo
!$acc end kernels

        !---
        !--- new formulation for cloud base and top
        !--- 1st guess for entraiment ( includes entrainment formulation ~ 1/z )
!$acc parallel
!$acc loop gang
        DO i=its,itf
            IF(ierr(I) /= 0) cycle
!$acc loop vector private(frh)
            do k=kts,ktf
                frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
                entr_rate_2d(i,k)=entr_rate*(1.3-frh)*max(min(1.,(qeso_cup(i,max(k,klcl(i)))&
            !                                                             /qeso_cup(i,klcl(i)))**3)   ,0.1) ! 
                                                                    /qeso_cup(i,klcl(i)))**1.25),0.1) ! P9_P6
            enddo
        enddo
!$acc end parallel

        !
        !--- Determine the level of convective cloud base (kbcon) and top (ktop),
        call cup_cloud_limits(cumulus,ierrc,ierr,cap_max_increment,cap_max,entr_rate &
                            ,heo_cup,heso_cup,qo_cup,qeso_cup,po,po_cup,zo_cup,heo,hkbo &
                            ,entr_rate_2d,hcot,k22,kbmax,klcl,kbcon,ktop            &
                            ,use_excess,zqexec,ztexec,xland &
                            ,itf,ktf,its,ite, kts,kte)

        !
        !--- Define detrainment
!$acc parallel loop gang
        DO i=its,itf
            if(ierr(i) /= 0)cycle
!$acc loop vector
            do k=kts,ktf
                cd(i,k)=0.75*entr_rate_2d(i,k)!+0.5e-3
            enddo
        enddo
!$acc end parallel

        !--- Get inversion layers for ktop
        ! Note : Current verification data does not test this branch
        IF(USE_INV_LAYERS) THEN
            call get_inversion_layers(cumulus,ierr,psur,po_cup,tn_cup,zo_cup,k_inv_layers,&
                                        dtempdz,itf,ktf,its,ite, kts,kte)
!$acc parallel loop vector
            do i=its,itf
                if(ierr(i) /= 0)cycle
                ktop(i) = min(ktop(i),k_inv_layers(i,shal)+1)
            enddo
!$acc end parallel
        ENDIF
        !
        !--- Check if ktop is above 700hPa layer for shallow convection
!$acc parallel 
!$acc loop vector
        do i=its,itf
            if(ierr(i) /= 0)cycle
            !print*,"sta=",Kbcon(i),kstabm(i),kstabi(i),p_cup(i,ktop(i)),z_cup(i,kstabi(i))
            if(po_cup(i,ktop(i)) < 650.) then
                ierr(i)=26
                ierrc(i)='shallow convection with cloud top above 650 hPa'
            endif
        enddo
!$acc loop vector
        do i=its,itf
            if(ierr(i) /= 0)cycle
            if(ktop(i) < kbcon(i)+2)then
                ierr(i)=5
                ierrc(i)='ktop too small'
            endif
        enddo
!$acc end parallel
        !
        !--- Determine the normalized mass flux for updraft
        !
        cumulus_1 = trim(cumulus)
        cumulus_2 = trim(cumulus)//"_up"
!$acc kernels
        do i=its,itf
            ! zuo(i,:)=0.
            if(ierr(i) /= 0) CYcle
            do k = kts,kte
                aa(k) = 0.0
                bb(k) = po_cup(i, k)
            enddo
            ! call get_zu_zd_pdf(trim(cumulus),trim(cumulus)//"_up",ierr(i),k22(i),ktop(i),zuo(i,kts:kte),kts,kte,ktf  &
            !                 ,kpbl(i),k22(i),kbcon(i),klcl(i),po_cup(i,kts:kte),psur(i),xland(i))
            call get_zu_zd_pdf( &
            cumulus_1,&
            cumulus_2, &
            ierr(i), &
            k22(i),&
            ktop(i),&
            ! zuo(i,kts:kte),&
            aa(kts:kte), &
            kts,&
            kte,&
            ktf,  &
            kpbl(i),&
            k22(i),&
            kbcon(i),&
            klcl(i),&
            ! po_cup(i,kts:kte),&
            bb(kts:kte), &
            psur(i),&
            xland(i))
            do k = kts,kte
                 zuo(i,k) = aa(k)
            enddo
        enddo
!$acc end kernels

!$acc parallel loop vector
        do i=its,itf
            if(ierr(i).eq.0)then
                !zuo(i,1)=0.
                xzu(i,1)= zuo(i,1)
                zu (i,1)= zuo(i,1)
                !~ do k=ktop(i)+1,ktf
                !~ zuo(i,k)=0.
                !~ enddo
            endif
        enddo
!$acc end parallel

        !---tmp srf
        !       do i=its,itf
        !         !IF(ierr(I).eq.0)THEN
        !          write(9,20)'ksh',ntimes,ierr(i),k22(i),kpbl(i),klcl(i),kbcon(i),ktop(i)
        !             20 format(1x,A5,7i4)
        !
        !         !endif
        !       enddo
        !---tmp srf
        !
        !--- calculate mass entrainment and detrainment
        !
        CALL get_lateral_massflux(itf,ktf, its,ite, kts,kte &
                                ,ierr,ktop,zo_cup,zuo,cd,entr_rate_2d        &
                                ,up_massentro, up_massdetro ,up_massentr, up_massdetr &
                                ,cumulus,kbcon,k22,up_massentru,up_massdetru,lambau)

!$acc kernels
        hc  =0.
        DBY =0.
        hco =0.
        DBYo=0.
        uc  =0.
        vc  =0.
        do i=its,itf
            IF(ierr(i).eq.0)THEN

                do k = kts,kte
                    aa(k) = po(i,k)
                    bb(k) = u_cup(i,k)
                    cc(k) = v_cup(i,k)
                enddo

                do k=kts,max(klcl(i)-1,kts)!kbcon(i)-1
                    hc (i,k)=hkb (i)
                    hco(i,k)=hkbo(i)
                    !-get uc and vc as average between layers below k22
                    call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),aa,bb,uc(i,k),k22(i))
                    call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),aa,cc,vc(i,k),k22(i))
                enddo
                k=klcl(i)!kbcon(i)
                hc  (i,k)=hkb (i)
                dby (i,k)=hkb (i)-hes_cup (i,k)
                hco (i,k)=hkbo(i)
                dbyo(i,k)=hkbo(i)-heso_cup(i,k)
                uc  (i,k)=uc  (i,max(k-1,kts))
                vc  (i,k)=vc  (i,max(k-1,kts))
            endif ! ierr
        enddo
!$acc end kernels

!$acc kernels
        do 42 i=its,itf
            if(ierr(i).eq.0)then
                do k=2,ktf-1
                    xzu(i,k)= zuo(i,k)
                    zu (i,k)= zuo(i,k)
                enddo
                do k=klcl(i)+1,ktop(i)
                    !do k=kbcon(i)+1,ktop(i)
                    hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1)+ &
                                                    up_massentr(i,k-1)*he(i,k-1))   /            &
                                (1.e-16+zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
                    dby(i,k)=hc(i,k)-hes_cup(i,k)

                    hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                                                        up_massentro(i,k-1)*heo(i,k-1))   /            &
                                (1.e-16+zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))

                    dbyo(i,k)=hco(i,k)-heso_cup(i,k)

                    denomU=1.e-16+(zu(i,k-1)-.5*up_massdetru(i,k-1)+up_massentru(i,k-1))
                    uc(i,k)=(uc(i,k-1)*zu(i,k-1)-.5*up_massdetru(i,k-1)*uc(i,k-1)+     &
                                                    up_massentru(i,k-1)*us(i,k-1)      &
                            -pgcon*.5*(zu(i,k)+zu(i,k-1))*(u_cup(i,k)-u_cup(i,k-1))) /  denomU

                    vc(i,k)=(vc(i,k-1)*zu(i,k-1)-.5*up_massdetru(i,k-1)*vc(i,k-1)+     &
                                                    up_massentru(i,k-1)*vs(i,k-1)      &
                            -pgcon*.5*(zu(i,k)+zu(i,k-1))*(v_cup(i,k)-v_cup(i,k-1))) /  denomU

                enddo
                !-
                do k=ktop(i)+1,ktf-1
                    ! HC(i,K)=hes_cup(i,k)
                    ! HCo(i,K)=heso_cup(i,k)
                    ! DBY(I,K)=0.
                    ! DBYo(I,K)=0.
                    zu (i,k)=0.
                    xzu(i,k)=0.
                    zuo(i,k)=0.
                    uc (i,k)= u_cup(i,k)
                    vc (i,k)= v_cup(i,k)
                enddo
                if(i.eq.ipr)then
                    write(0,*)'hcnew = '
                    do k=kts,ktf
                        write(0,*)k,hco(i,k),dbyo(i,k)
                    enddo
                endif
            endif
        42    continue
!$acc end kernels

        call cup_up_moisture(cumulus,klcl,ierr,zo_cup,qco,qrco,pwo,pwavo,hco,tempco,xland, &
            po,p_cup,kbcon,ktop,cd,dbyo,clw_all, &
            t_cup,qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,        &
            ZQEXEC,use_excess,ccn,rho,up_massentr,up_massdetr,psum,psumh,c1d,&
            1,itf,ktf,ipr,jpr,its,ite, kts,kte)

!$acc kernels
        do i=its,itf
            if(ierr(i).eq.0)then
                cupclw(i,kts:ktf)=qrco(i,kts:ktf)
                qco (i,ktop(i)+1:ktf-1)=QESo_cup(I,ktop(i)+1:ktf-1)
            endif
        enddo
!$acc end kernels

        !--- calculate in-cloud air temperature for CAPE  and vertical velocity
        !
!$acc kernels
        do i=its,itf
            if(ierr(i) == 0)then
                do k=kts,ktf
                    tempco(i,k) = (1./cp)*(hco(i,k)-g*zo_cup(i,k)-xlv*qco(i,k))
                    !~ print*,"tempco",k,tempco(i,k),t_cup(i,k),zo_cup(i,k)
                enddo
            else
                tempco(i,:)=tn_cup(i,:)
            endif
        enddo
!$acc end kernels

        !
        !--- calculate vertical velocity
        !
        call cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd,zo,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup &
                        ,tempco,qco,qrco,qo &
                        ,klcl,kbcon,ktop             &
                        ,ierr,itf,ktf,its,ite, kts,kte)
        !
        !--- calculate cape for updrafts
        !
        call cup_up_cape(cape,z,zu,dby,GAMMA_CUP,t_cup, &
                        k22,kbcon,ktop,ierr,tempco,qco,qrco, qo_cup  , &
                        itf,ktf,its,ite, kts,kte)

        !
        !
        !--- get the environmental mass flux
        !
        ! do i=its,itf
        !     zenv(i,:) = 0.0
        !     if(ierr(i) /= 0) cycle
        !     zenv(i,:) = zuo(i,:)
        ! enddo
!$acc parallel loop gang vector collapse(2)
        do k = kts,kte
            do i = its,ite
                zenv(i,k) = 0.0
                if(ierr(i).eq.0) zenv(i,k) = zuo(i,k)
            enddo
        enddo
!$acc end parallel

        !
        !--- change per unit mass that a model cloud would modify the environment
        !
        !--- 1. in bottom layer
        !
!$acc parallel loop gang vector collapse(2)
        do k=kts,kte
            do i=its,itf
                dellah (i,k)=0.
                dellaq (i,k)=0.
                dellaqc(i,k)=0.
                dellu  (i,k)=0.
                dellv  (i,k)=0.
            enddo
        enddo
!$acc end parallel
!
!----------------------------------------------  cloud level ktop
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level k+2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
!
!----------------------------------------------  cloud level k+1
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k
!
!----------------------------------------------  cloud level k
!
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!-------------------------------------------z_cup(3) cloud level 3
!
!- - - - - - - - - - - - - - - - - - - - - -z    (2) model level 2
!
!-------------------------------------------z_cup(2) cloud level 2
!
!- - - - - - - - - - - - - - - - - - - - - -z    (1) model level 1
!
!-------------------------------------------z_cup(1) cloud level 1-SURFACE LEVEL
        IF(VERT_DISCR == 1) THEN
!$acc parallel loop gang
            do i=its,itf
                if(ierr(i).eq.0)then
                    !trash=0.
                    !trash2=0.
!$acc loop vector private(entup, detup, totmas, dp, C_up)
                    do k=kts,ktop(i)

                        ! entrainment/detrainment for updraft
                        entup=up_massentro(i,k)
                        detup=up_massdetro(i,k)

                        totmas=detup-entup+zuo(i,k+1)-zuo(i,k)
                        if(abs(totmas).gt.1.e-6)then
                            print*,'BLAH1!'
                            ! write(0,*)'****shallow i k error',i,k,totmas
                            ! write(0,*)'k22-kbcon-ktop',k22(i),kbcon(i),ktop(i)
                            ! write(0,123)detup,entup,zuo(i,k+1),zuo(i,k)
                            ! 123     format(1X,i2,5E12.4)
                        endif
                        dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                        dellu (i,k) =-(zuo(i,k+1)*(uc (i,k+1)-u_cup(i,k+1) ) -      &
                                        zuo(i,k  )*(uc (i,k  )-u_cup(i,k  ) ) )*g/dp

                        dellv (i,k) =-(zuo(i,k+1)*(vc (i,k+1)-v_cup(i,k+1) ) -      &
                                        zuo(i,k  )*(vc (i,k  )-v_cup(i,k  ) ) )*g/dp

                        dellah(i,k) =-(zuo(i,k+1)*(hco(i,k+1)-heo_cup(i,k+1) )-     &
                                        zuo(i,k  )*(hco(i,k  )-heo_cup(i,k  ) ))*g/dp

                        !-- take out cloud liquid water for detrainment
                        dellaqc(i,k)=   detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

                        !-- condensation source term = detrained + flux divergence of
                        !-- cloud liquid water (qrco)
                        C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                                zuo(i,k  )* qrco(i,k  )  )*g/dp

                        !-- water vapor budget (flux divergence of Q_up-Q_env - condensation term)
                        dellaq(i,k) =-(zuo(i,k+1)*(qco(i,k+1)-qo_cup(i,k+1) ) -         &
                                        zuo(i,k  )*(qco(i,k  )-qo_cup(i,k  ) ) )*g/dp &
                                        - C_up

                        !- check water conservation liq+condensed
                        !trash =trash + (dellaq(i,k)+dellaqc(i,k))*dp/g
                        !trash2=trash2+  dellah(i,k)*dp/g

                    enddo   ! k

                    !~ write(0,*)'=>checking Q/H cons for shallow= ',trash,trash2
                    !~ if(abs(trash)>1.e-6 .or. abs(trash2)>1.e-6) then
                    !~ write(6,*)'=> not Q/H cons for shallow= ',i,trash,trash2
                    !~ stop "shallow not conserving"
                    !~ endif

                endif
            enddo
!$acc end parallel
        ELSEIF(VERT_DISCR == 2) THEN
!$acc parallel loop gang
            do i=its,itf
                if(ierr(i) /= 0) CYCLE
                !trash=0.
                !trash2=0.
!$acc loop vector private(entup, detup, totmas, dp, C_up)
                do k=kts,ktop(i)

                    ! entrainment/detrainment for updraft
                    entup=up_massentro(i,k)
                    detup=up_massdetro(i,k)

                    totmas=detup-entup+zuo(i,k+1)-zuo(i,k)
                    if(abs(totmas).gt.1.e-6)then
                        print*, 'BLAH2!'
                        ! write(0,*)'****shallow i k error',i,k,totmas
                        ! write(0,*)'k22-kbcon-ktop',k22(i),kbcon(i),ktop(i)
                        ! write(0,125)detup,entup,zuo(i,k+1),zuo(i,k)
                        ! 125     format(1X,i2,5E12.4)
                    endif
                    dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                    dellu (i,k) =-(zuo(i,k+1)*(uc (i,k+1)-us(i,k+1) ) -      &
                                    zuo(i,k  )*(uc (i,k  )-us(i,k  ) ))*g/dp

                    dellv (i,k) =-(zuo(i,k+1)*(vc (i,k+1)-vs(i,k+1) ) -      &
                                    zuo(i,k  )*(vc (i,k  )-vs(i,k  ) ))*g/dp

                    !------------
                    !           dellah(i,k) =-(zuo(i,k+1)*(hco(i,k+1)-heo(i,k+1) )-      &
                    !                          zuo(i,k  )*(hco(i,k  )-heo(i,k  ) ))*g/dp
                    !--using vert_discr=1 for H
                        dellah(i,k) =-(zuo(i,k+1)*(hco(i,k+1)-heo_cup(i,k+1) )-     &
                                        zuo(i,k  )*(hco(i,k  )-heo_cup(i,k  ) ))*g/dp
                    !-------------

                    !-- take out cloud liquid water for detrainment
                    dellaqc(i,k)=   detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

                    !-- condensation source term = detrained + flux divergence of
                    !-- cloud liquid water (qrco)
                    C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                            zuo(i,k  )* qrco(i,k  )  )*g/dp

                    !-- water vapor budget (flux divergence of Q_up-Q_env - condensation term)
                    dellaq(i,k) =-(zuo(i,k+1)*(qco(i,k+1)-qo(i,k+1) ) -         &
                                    zuo(i,k  )*(qco(i,k  )-qo(i,k  ) ) )*g/dp    &
                                    - C_up

                    !- check water conservation liq+condensed
                    !trash =trash + (dellaq(i,k)+dellaqc(i,k))*dp/g
                    !trash2=trash2+  dellah(i,k)*dp/g

                enddo   ! k

                !write(0,*)'=>checking Q/H cons for shallow= ',trash,trash2
                !if(abs(trash)>1.e-6 .or. abs(trash2)>1.e-6) then
                !  write(6,*)'=> not Q/H cons for shallow= ',i,trash,trash2
                !  stop "shallow not conserving"
                !endif
            enddo
!$acc end parallel
        ENDIF
        !
        !--- using dellas, calculate changed environmental profiles
        !
        mbdt=3.e-4
        !do k=kts,ktf
        !  write(0,*)'zuo,k22 = ',k,zuo(1,k),k22(1)
        !enddo

!$acc parallel loop gang vector collapse(2)
        do k=kts,ktf
            do i=its,itf
                dellat(i,k)=0.
                if(ierr(i).eq.0)then
                    !
                    XHE(I,K)=(DELLAH(I,K))*MBDT+HEO(I,K)
                    XQ (I,K)=(DELLAQ(I,K)+DELLAQC(i,k))*MBDT+QO(I,K)
                    XQ (I,K)=max(XQ(I,K),1.E-08)
                    !- do not feed dellat with dellaqc if
                    !- the detrainment of liquid water will be used as
                    !- a source for cloud microphysics (then microphysics will
                    !- evaporate this excess and provide the cooling at the
                    !- detrainment region)
                    if(COUPL_MPHYSICS) then
                        DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xlv*(DELLAQ(I,K)))
                    else
                        DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xlv*(DELLAQ(I,K)+DELLAQC(i,k)))
                        !-adding dellaqc to dellaq:
                        DELLAQ (I,K)= DELLAQ(I,K)+DELLAQC(I,K)
                        DELLAQC(I,K)= 0.0
                    endif
                    XT(I,K)=((1./cp)*DELLAH(i,k)-(xlv/cp)*(DELLAQ(i,k)+ DELLAQC(i,k)))*MBDT+TN(I,K)
                ENDIF
            enddo
        enddo
!$acc end parallel

        !
        ! now for shallow forcing
        !
!$acc parallel loop gang private(xff_shal, blqe, trash, tcold, dp, thot, efic, fin, fsum, fixouts)
        do i=its,itf
            xmb(i)=0.
            xff_shal(1:9)=0.
            if(ierr(i).eq.0)then

                !         xmbmax(i)=0.1
                xmbmax(i)=100.*(po(i,kbcon(i))-po(i,kbcon(i)+1))/(g*dtime)
                !
                !- limiting the mass flux at cloud base
                xmbmax(i)=min(xmbmaxshal,xmbmax(i))

                !
                !- closure from Grant (2001)
                xff_shal(1)=.03*zws(i)*rho(i,kpbl(i))
                xff_shal(2)=.03*zws(i)*rho(i,kpbl(i))
                xff_shal(3)=.03*zws(i)*rho(i,kpbl(i))
                !----
                !- closure from boundary layer QE (Raymond 1995)
                blqe=0.
                trash=0.
                if(k22(i).lt.kpbl(i)+1)then
                    do k=kts,kbcon(i)
                    blqe=blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
                    enddo
                    !trash=max((hc(i,kbcon(i))-he_cup(i,kbcon(i))),1.e1)
                    trash = max((hco(i,kbcon(i))-heo_cup(i,kbcon(i))),1.e1)
                    xff_shal(7)=max(0.,blqe/trash)
                    !print*,"blqe=", xff_shal(7),blqe,trash
                else
                    xff_shal(7)=0.0
                endif
                xff_shal(8)= xff_shal(7)
                xff_shal(9)= xff_shal(7)
                !----
                !- closure from the heat-engine principle
                !- Renno and Ingersoll(1996), Souza et al (1999)
                    !- get the averaged environment temperature between cloud base
                    !- and cloud top
                tcold=0.
                do k=kbcon(i),ktop(i)
                    dp   = po_cup(i,k)-po_cup(i,k+1)
                    tcold= tcold + t_cup(i,k)*dp
                enddo
                tcold=tcold/(po_cup(i,kbcon(i))-po_cup(i,ktop(i)+1))

                !-surface temperature
                thot=tsur(i)  ! + ztexec(i)
                !- thermodynamic eficiency
                efic = max(0.05, (thot-tcold)/thot )
                !efic = max(0.0, (thot-tcold)/thot )

                !- total heat flux from surface
                fin = max(0.0, h_sfc_flux(i)+le_sfc_flux(i))

                !- mass flux at cloud base
                !if(cape(i) > 0.0 .and. h_sfc_flux(i) >0.0 ) then
                if(cape(i) > 0.0  ) then
                    xff_shal(4) = efic * fin / cape(i)
                    !if(xff_shal(4)>0.0001) then
                    !~ print*,"xff_shal1=",i,h_sfc_flux(i),le_sfc_flux(i),tsur(i)
                    !~ print*,"xff_shal2=",xff_shal(4) , efic , fin , cape(i)
                    !~ flush(6)
                    !endif
                else
                    xff_shal(4) = 0.0
                endif
                xff_shal(5)=xff_shal(4)
                xff_shal(6)=xff_shal(4)

                fsum=0.
                do k=1,9
                    !- heat engine closure is providing too low values for mass fluxes.
                    !- until this is checked, the ensemble closure will be calculatd
                    !- only using the closures BLQE and Wstar
                    if(k.ge.4 .and. k.le.6) cycle
                    xmb(i)=xmb(i)+xff_shal(k)
                    fsum=fsum+1.
                enddo
                !- ensemble average of mass flux
                xmb(i)=min(xmbmax(i),xmb(i)/fsum)

                if(ichoice > 0) xmb(i)=min(xmbmax(i),xff_shal(ichoice))

                !~ print*,'err,xffs',ierr(i),xff_shal(1:9),xmb(i),xmbmax(i)
                if(xmb(i).eq.0.)then
                    ierr(i)=22
                    ierrc(i)='22'
                endif
                if(xmb(i).lt.0.)then
                    ierr(i)=21
                    ierrc(i)='21'
                    !write(0,*)'neg xmb,i,xmb for shallow = ',i,k22(i),ierr(i)
                endif
            endif

            !~ ! max allowable tendency over tendency that would lead to too small mix ratios
            !~ !
            !~ trash=(1.e-12 -q(i,k)) /((dsubq(i,k)+dellaq(i,k))*dtime)
            !~ xmb(i)=(1.e-12 -q(i,k))/((dsubq(i,k)+dellaq(i,k))*dtime)


            !--- check outtem and outq for high/unphysical values
            !--- criteria: if abs (dT/dt or dQ/dt) > 100 K/day => fix xmb
            IF(ierr(i) == 0) then
                fixouts=xmb(i) *86400.*max(maxval(abs(dellat(i,kts:ktop(i)))),&
                                    (xlv/cp)*maxval(abs(dellaq(i,kts:ktop(i)))) )

                if(fixouts > 100.) then ! K/day
                    fixouts=100./(fixouts)
                    xmb(i) = xmb(i)*fixouts
                endif
                !
                ! final tendencies
                !
                do k=kts,ktop(i)
                    outt (i,k)=dellat (i,k)*xmb(i)
                    outq (i,k)=dellaq (i,k)*xmb(i)
                    outqc(i,k)=dellaqc(i,k)*xmb(i)
                    outu (i,k)= dellu (i,k)*xmb(i)
                    outv (i,k)= dellv (i,k)*xmb(i)
                enddo
            ENDIF
        enddo
!$acc end parallel

!----------------------------------------------------------------
!--- cloud fraction calculation - not in use --------------------
!  do i=its,itf
!    clfrac(:,:)=0.
!    if(ierr(i) == 0) then
!     do k=kts,ktop(i)
!       frh = min(q_cup(i,k)/qes_cup(i,k),0.999)
!       clfrac(i,k)=(frh**0.25)*(1.-exp(-100.*qrco(i,k)/&
!                                ((1.-frh)*qes_cup(i,k))**0.49))
!       print*,"k=",clfrac(i,k),frh,qrco(i,k)
!     enddo
!   endif
!  enddo
!  do i=its,itf
!      !clfrac(i,:)=0.
!      if(ierr(i) /= 0) cycle
!       dummy1(kts:ktf) = xmb(i)* zuo(i,kts:ktf)
!       dummy2(kts:ktf) = 100.*po_cup(i,kts:ktf)
!       call get_cloud_fraction(ktf,kts,ktf                                              &
!                 ,dummy2(kts:ktf),zo_cup(i,kts:ktf),tn_cup(i,kts:ktf),qo_cup(i,kts:ktf) &
!                 ,qco (i,kts:ktf),  qrco(i,kts:ktf),  dummy1(kts:ktf),clfrac(i,kts:ktf) )
!  enddo
!----------------------------------------------------------------
!$acc parallel loop gang
        do i=its,itf
            if(ierr(i) == 0) then
!$acc loop vector
                do k=kts,ktop(i)
                    !clwup(i,k) = qrco(i,k) ! ice/liquid water
                    tup   (i,k) = tempco(i,k) !updraft temp
                enddo
            endif
        enddo
!$acc end parallel

        !
        !--------------------------------------------------------------------------------------------!
        !--------------------------------------------------------------------------------------------!
        !- section for atmospheric composition
        !--------------------------------------------------------------------------------------------!
        IF(USE_TRACER_TRANSP==1) THEN

            !-- for debug only
            !if(jl==1) then 
            !  se_chem_update(1,:,:) = se_chem(1,:,:)
            !else
            !  se_chem(1,:,:) = se_chem_update(1,:,:)
            !endif   
            !do i=its,itf
            !     if(ierr(i) /= 0) cycle
            !      massi=0.
            !      do k=kts,ktf!ktop(i)
            !            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            !            massi  =   massi + se_chem(1,i,k)*dp/g
            !      enddo
            !enddo
            !-- for debug only

            !-0) convert mass fluxes

            ! do i=its,itf
            !     if(ierr(i) /= 0) cycle
            !     zuo         (i,:) = xmb(i)* zuo(i,:)
            !     up_massentro(i,:) = xmb(i)* up_massentro(i,:)
            !     up_massdetro(i,:) = xmb(i)* up_massdetro(i,:)
            !     zenv        (i,:) = xmb(i)* zenv        (i,:)
            ! enddo

!$acc parallel loop gang vector
            do k = kts,kte
                do i = its,ite
                    if(ierr(i) == 0) then
                        zuo         (i,k) = xmb(i)* zuo(i,k)
                        up_massentro(i,k) = xmb(i)* up_massentro(i,k)
                        up_massdetro(i,k) = xmb(i)* up_massdetro(i,k)
                        zenv        (i,k) = xmb(i)* zenv        (i,k)
                    endif
                enddo
            enddo
!$acc end parallel


            !- 1) get mass mixing ratios at the cloud levels

            call cup_env_clev_chem(mtp,se_chem,se_cup_chem,ierr,itf,ktf,its,ite, kts,kte)

            !-2) determine the boundary condition for in-cloud tracer mixing ratios
            !
            !  do i=its,itf
            !     if(ierr(i) /= 0)cycle
            !     x_add = 0.
            !     do ispc=1,mtp
            !      call get_cloud_bc(kts,kte,po(i,kts:kte),se_cup_chem(ispc,i,kts:kte),sc_up_chem_b(ispc,i),k22(i),x_add)
            !     enddo
            !  enddo
            !-3) determine the in-cloud tracer mixing ratios
            !
!$acc kernels
            sc_up_chem=se_cup_chem
!$acc end kernels

!$acc parallel loop gang
            do i=its,itf
                if(ierr(i) /= 0) cycle
                start_level(i) = k22(i)
                !start_level(i) = kts
!$acc loop vector
                do k=kts,start_level(i)
                    !sc_up_chem(:,i,k)=sc_up_chem_b(:,i)
                    sc_up_chem(:,i,k)=se_cup_chem(:,i,k)
                enddo
            enddo
!$acc end parallel

!$acc parallel loop gang
            do i=its,itf
                if(ierr(i) /= 0) cycle
!$acc loop seq
                do k=start_level(i)+1,ktop(i)+1
                    !-- entr,detr, mass flux ...
                    XZZ=             zuo(i,k-1)
                    XZD=0.5*up_massdetro(i,k-1)
                    XZE=    up_massentro(i,k-1)
                    denom =  (XZZ-XZD+XZE)
                    if(denom > 0.) then
                        !-- transport + mixing
                        sc_up_chem(:,i,k) = (sc_up_chem(:,i,k-1)*XZZ - sc_up_chem(:,i,k-1)*XZD + XZE*se_chem(:,i,k-1)) / denom
                    else
                        sc_up_chem(:,i,k) =  sc_up_chem(:,i,k-1)
                    endif
                enddo
            enddo
!$acc end parallel
            !
            !-4) determine the vertical transport 
            !
            !---a) change per unit mass that a model cloud would modify the environment
!$acc kernels
            do i=its,itf
                if(ierr(i) /= 0) cycle
    
                ! if(USE_FLUX_FORM == 1 .and. alp1 == 0. ) then !- flux form + source/sink terms + time explicit + FCT

                !     if(use_fct == 0 ) then 
                !         do k=kts,ktop(i)
                !             dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                !             out_chem(:,i,k) =-(zuo(i,k+1)*(sc_up_chem(:,i,k+1)-se_cup_chem(:,i,k+1) ) -                 &
                !                                 zuo(i,k  )*(sc_up_chem(:,i,k  )-se_cup_chem(:,i,k  ) ))*g/dp             
                !         enddo
        
                !     else 

                !         !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
                !         sub_tend  = 0.
                !         trcflx_in = 0. 
                !         dtime_max = dtime
                !         massflx  (i,kts)=0.

                !         do k=kts+1,ktop(i)+1
                !             dp  	       = 100.*(po_cup(i,k)-po_cup(i,k+1))
                !             trcflx_in (:,k) =-zuo(i,k) *se_cup_chem(:,i,k)  !* xmb(i)
                !             massflx   (i,k) =-zuo(i,k) 		       !* xmb(i)
                !             dtime_max=min(dtime_max,.5*dp) 
                !         enddo
                !         !--if dtime_max<dtime => needs a loop to update from t to t+dtime (check this!)
                !         !if( dtime_max < dtime ) stop "dtime_max < dtime in GF scheme"
                        
                !         do ispc=1,mtp
                !             call fct1d3 (ktop(i),kte,dtime_max,po_cup(i,:),se_chem(ispc,i,:),massflx(i,:),trcflx_in(ispc,:),sub_tend(ispc,:))
                !         enddo

                !         do k=kts,ktop(i)
                !             dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                !             out_chem(:,i,k) =-(zuo(i,k+1)*(sc_up_chem(:,i,k+1)) - zuo(i,k)*(sc_up_chem(:,i,k)))*g/dp  	

                !             !- update with subsidence term from FCT scheme
                !             out_chem(:,i,k) = out_chem(:,i,k) + sub_tend(:,k)
                !         enddo
                !     endif
                ! endif
   

                if(USE_FLUX_FORM == 1 .and. alp1 > 0. ) then !- flux form + source/sink terms + time explicit/implicit + upstream
                    alp0=1.-alp1

                    do k=kts,ktop(i)
                        dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                        beta1 = dtime*g/dp
                        bb(k) = 1.+alp1*beta1*zenv(i,k)
                        cc(k) =   -alp1*beta1*zenv(i,k+1)
                        
                        ddtr(:,k) = se_chem(:,i,k) - (zuo(i,k+1)*sc_up_chem(:,i,k+1) - zuo(i,k)*sc_up_chem(:,i,k))*beta1 
                            
                        ddtr(:,k) = ddtr(:,k) + alp0*beta1*(-zenv(i,k)*se_chem(:,i,k) +zenv(i,k+1)*se_chem(:,i,k+1))
                    
                    enddo
                    do ispc = 1, mtp
                        do k = kts,ktop(i)
                            ddu(k) = ddtr(ispc, k)
                        enddo
                        ! call bidiag (ktop(i),bb (kts:ktop(i)),cc (kts:ktop(i)),ddtr (ispc,kts:ktop(i)))
                        call bidiag (ktop(i),bb (kts:ktop(i)),cc (kts:ktop(i)),ddu(kts:ktop(i)))
                        do k = kts,ktop(i)
                            ddtr(ispc, k) = ddu(k)
                        enddo
                        out_chem(ispc,i,kts:ktop(i))=(ddtr(ispc,kts:ktop(i))-se_chem(ispc,i,kts:ktop(i)))/dtime
                    enddo
    
                endif
  

            enddo
!$acc end kernels
            !-- for debug only
            !   do i=its,itf
            !        massf = 0.
            !        if(ierr(i) /= 0) cycle
            !        do k=kts,ktf !ktop(i)+1
            !            se_chem_update(1,i,k) = se_chem_update(1,i,k) + out_chem(1,i,k) * dtime
            !            se_chem_update(1,i,k) = max(0.,se_chem_update(1,i,k))
            !            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            !            massf = massf + se_chem_update(1,i,k)*dp/g
            !        enddo
            !        !if(abs((massf-massi)/(1.e-12+massi))>1.e-6) then 
            !          print*,"mass con=>",(massf-massi)/(1.e-12+massi)
            !          !stop "MASS CON"
            !        !endif
            !   enddo
            !19 FORMAT(1x,I3,1x,5E14.3);18 FORMAT(1x,I3,1x,4E14.3);20 FORMAT(1x,I3,1x,11E16.6)
            !-- for debug only
            !-4) convert back mass fluxes

            ! do i=its,itf
            !     if(ierr(i) /= 0) cycle
            !     zuo         (i,:) = zuo(i,:)          / (xmb(i) + 1.e-16)                                          
            !     up_massentro(i,:) = up_massentro(i,:) / (xmb(i) + 1.e-16)        
            !     up_massdetro(i,:) = up_massdetro(i,:) / (xmb(i) + 1.e-16)         
            !     zenv        (i,:) = zenv        (i,:) / (xmb(i) + 1.e-16)        
            ! enddo

!$acc parallel loop gang vector
            do k = kts,kte
                do i = its,ite
                    if(ierr(i) == 0) then
                        zuo         (i,k) = zuo(i,k)          / (xmb(i) + 1.e-16)                                          
                        up_massentro(i,k) = up_massentro(i,k) / (xmb(i) + 1.e-16)        
                        up_massdetro(i,k) = up_massdetro(i,k) / (xmb(i) + 1.e-16)         
                        zenv        (i,k) = zenv        (i,k) / (xmb(i) + 1.e-16) 
                    endif
                enddo
            enddo
!$acc end parallel
        !
        !--------------------------------------------------------------------------------------------!
        ENDIF !- end of section for atmospheric composition
        !--------------------------------------------------------------------------------------------!
!$acc end data
        !- begin: for GATE soundings
        if(wrtgrads) then
            cty='2'
            ! do i=its,itf
            !     do k=kts,kte!max(1,ktop(i))
            !         nvar=101
            !         call set_grads_var(jl,k,nvar,out_chem(1,i,k)*86400,"outchem"//cty ,' outchem','3d')
            !         call set_grads_var(jl,k,nvar,sc_up_chem(1,i,k),"sc_chem"//cty ,' sc_up_chem','3d')
            !         call set_grads_var(jl,k,nvar,se_chem(1,i,k),"se_chem"//cty ,' se_chem','3d')

            !         call set_grads_var(jl,k,nvar,zuo(i,k),"zup"//cty,'norm m flux up sh','3d')
            !         call set_grads_var(jl,k,nvar,xmb(i),"xmb"//cty,'base m flux up sh (kg/s/m^2)','2d')
            !         call set_grads_var(jl,k,nvar,dellah(i,k),"delh"//cty ,'dellah_sh','3d')
            !         call set_grads_var(jl,k,nvar,dellaq(i,k)*86400.*xlv/cp*xmb(i), "dellq"//cty ,'dellaq_sh','3d')
            !         call set_grads_var(jl,k,nvar,dellaqc(i,k)*86400.*xlv/cp*xmb(i),"dellqc"//cty ,'dellaqc_sh','3d')
            !         call set_grads_var(jl,k,nvar,up_massentro(i,k),"upent"//cty ,'up_massentr_sh','3d')
            !         call set_grads_var(jl,k,nvar,up_massdetro(i,k),"updet"//cty ,'up_massdetr_sh','3d')
            !         call set_grads_var(jl,k,nvar,outt(i,k)*86400.,"outt"//cty ,'outt_sh K/s','3d')
            !         call set_grads_var(jl,k,nvar,outu(i,k)*86400.,"outu"//cty ,'outu_sh m/s/s','3d')
            !         call set_grads_var(jl,k,nvar,outv(i,k)*86400.,"outv"//cty ,'outv_sh m/s/s','3d')
            !         call set_grads_var(jl,k,nvar,outq(i,k)*86400.*xlv/cp,"outq"//cty ,'outq_sh K/s','3d')
            !         call set_grads_var(jl,k,nvar,outqc(i,k)*86400.*xlv/cp,"outqc"//cty ,'outqc_sh K/s','3d')
            !         call set_grads_var(jl,k,nvar,float(ierr(i)),"ierr"//cty ,'ierrsh #','2d')
            !         !nvar=70
            !         !frh = (qeso_cup(i,k)/qeso_cup(i, max(1,kbcon(i))) )
            !         ! frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
            !         ! frh= entr_rate*(1.3-frh)*min(1.,(qeso_cup(i,k)/qeso_cup(i,max(1,kbcon(i))))*0.125)
            !         !call set_grads_var(jl,k,nvar,frh,"frh"//cty ,'frh #','3d')

            !         call set_grads_var(jl,k,nvar,entr_rate_2d(i,k),"entr"//cty ,' m-1','3d')
            !         call set_grads_var(jl,k,nvar,cd(i,k),"detr"//cty ,' m-1','3d')
            !         !~ call set_grads_var(jl,k,nvar,DER_Q1(i,k),"dq1"//cty ,'dq1 zu*q*g/dp','3d')
            !         !~ call set_grads_var(jl,k,nvar,DER_Q2(i,k),"dq2"//cty ,'dq2 zu*q*g/dp','3d')
            !         !~ call set_grads_var(jl,k,nvar,QCUP(i,k),"qcup"//cty ,'C_UP','3d')
            !         !~ call set_grads_var(jl,k,nvar,QEDN(i,k),"qedn"//cty ,'E_DN','3d')
            !         !~ call set_grads_var(jl,k,nvar,DIF_Q1(i,k),"df1"//cty ,'df1 qco-qo','3d')
            !         !~ call set_grads_var(jl,k,nvar,DIF_Q2(i,k),"df2"//cty ,'df2 qco-qo','3d')
            !         !nvar=80

            !         call set_grads_var(jl,k,nvar,t_cup(i,k)-273.15,"te"//cty ,' K','3d')
            !         call set_grads_var(jl,k,nvar,1000.*q_cup(i,k),"qe"//cty ,' kg kg-1','3d')
            !         call set_grads_var(jl,k,nvar,heo_cup(i,k),"heo"//cty ,' he','3d')
            !         call set_grads_var(jl,k,nvar,vvel2d(i,k),"W2d"//cty ,' vvel','3d')
            !         call set_grads_var(jl,k,nvar,vvel1d(i),"W1d"//cty ,' vvel','2d')
            !         call set_grads_var(jl,k,nvar,he_cup(i,k),"he"//cty ,' he','3d')
            !         call set_grads_var(jl,k,nvar,hcot(i,k),"hcot"//cty ,' hcot','3d')
            !         call set_grads_var(jl,k,nvar,heso_cup(i,k),"heso"//cty ,' he','3d')
            !         call set_grads_var(jl,k,nvar,qeso_cup(i,k),"qeso"//cty ,' qes','3d')
            !         call set_grads_var(jl,k,nvar,hes_cup(i,k),"hes"//cty ,' he','3d')
            !         call set_grads_var(jl,k,nvar,dbyo(i,k),"buoy"//cty ,' buoy','3d')
            !         call set_grads_var(jl,k,nvar,HKBo(i),"hkbo"//cty ,' H','2d')
            !         call set_grads_var(jl,k,nvar,HKB(i),"hkb"//cty ,' H','2d')
            !         call set_grads_var(jl,k,nvar,float(k22(i)),"k22"//cty ,' k22','2d')
            !         call set_grads_var(jl,k,nvar,z_cup(i,kbcon(i)),"zbcon"//cty ,' m','2d')
            !         call set_grads_var(jl,k,nvar,z_cup(i,ktop(i)),"ztop"//cty ,' m','2d')
            !         call set_grads_var(jl,k,nvar,z_cup(i,klcl(i)),"zlcl"//cty ,' m','2d')
            !         call set_grads_var(jl,k,nvar,z_cup(i,kpbl(i)),"zpbl"//cty ,' m','2d')
            !         call set_grads_var(jl,k,nvar,float(kbcon(i)),"kbcon"//cty ,' m','2d')
            !         call set_grads_var(jl,k,nvar,float(ktop(i)),"ktop"//cty ,' m','2d')
            !         call set_grads_var(jl,k,nvar,float(klcl(i)),"klcl"//cty ,' m','2d')
            !         call set_grads_var(jl,k,nvar,float(kpbl(i)),"kpbl"//cty ,' m','2d')
            !         call set_grads_var(jl,k,nvar,zws(i),"ws"//cty ,' m/s','2d')
            !         !~ call set_grads_var(jl,k,nvar,1000.*dtempdz(i,k),"dtdz"//cty ,' K km-1','3d')
            !         !~ call set_grads_var(jl,k,nvar,clfrac(i,k),"clfrac"//cty ,'shcf #','3d')
            !     enddo
            !     call wrt_bin_ctl(1,kte,po(1,1:kte),cumulus)
            ! enddo
        endif
    !- end  : for GATE soundings-------------------------------------------

    ! done shallow
    END SUBROUTINE CUP_gf_sh

    SUBROUTINE CUP_gf_GEOS5(its,ite,kts,kte ,itf,ktf, mtp &
        ,cumulus           &
        ,ichoice           &
        ,entr_rate_input   &
        ,use_excess        &
        !input data
        ,dx                &
        ,stochastic_sig    &
        ,dtime             &
        ,kpbl              &
        ,ztexec            &
        ,zqexec            &
        ,ccn               &
        ,rho               &
        ,omeg              &
        ,t                 &
        ,q                 &
        ,z1                &
        , h_sfc_flux       &
        ,le_sfc_flux       &
        ,xlons             &
        ,xlats             &
        ,xland             &
        ,tn                &
        ,qo                &
        ,tn_bl             &
        ,qo_bl             &
        ,zo                &
        ,po                &
        ,tsur              &
        ,psur              &
        ,us                &
        ,vs                &
        ,dm2d              &
        ,se_chem           &
        ,zws               &
        ,dhdt              &
        ,last_ierr         &
        !output data
        ,outt              &
        ,outq              &
        ,outqc             &
        ,outu              &
        ,outv              &
        ,out_chem          &
        !- for convective transport
        ,ierr              &
        ,jmin              &
        ,klcl              &
        ,k22               &
        ,kbcon             &
        ,ktop              &
        ,kstabi            &
        ,kstabm            &
        ,pre               &
        ,xmb               &
        ,edto              &
        ,pwavo             &
        ,sig               &
        ,po_cup            &
        ,up_massentro      &
        ,up_massdetro      &
        ,dd_massentro      &
        ,dd_massdetro      &
        ,zuo               &
        ,zdo               &
        ,pwo               &
        ,pwdo              &
        ,qrco              &
        ,tup               &
        ,clfrac            &
        !- for convective transport-end
        !- for debug/diag
        ,AA0_,AA1_,AA2_,AA3_,AA1_BL_,AA1_CIN_,TAU_BL_,TAU_EC_   &
        ,revsu_gf          &
        ,prfil_gf           &
        )
        IMPLICIT NONE

        !-local settings
        LOGICAL, PARAMETER:: USE_LCL       =.FALSE.
        LOGICAL, PARAMETER:: USE_INV_LAYERS=.TRUE.

        CHARACTER*(*),INTENT(IN) :: cumulus
        INTEGER      ,INTENT(IN) :: itf,ktf,its,ite,kts,kte,ichoice,use_excess,mtp
        INTEGER      ,INTENT(IN),  DIMENSION (its:ite) ::   kpbl,last_ierr
        !
        ! outtem = output temp tendency (per s)
        ! outq   = output q tendency (per s)
        ! outqc  = output qc tendency (per s)
        ! pre    = output precip
        REAL,    DIMENSION (its:ite,kts:kte) ,INTENT (INOUT)   ::       &
            outu,outv,outt,outq,outqc,revsu_gf,prfil_gf

        REAL,    DIMENSION (its:ite)         ,INTENT (OUT  )   ::       &
            pre,sig

        !
        ! basic environmental input includes
        ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
        ! convection for this call only and at that particular gridpoint
        !
        REAL,    DIMENSION (its:ite,kts:kte)       ,INTENT (IN   )    ::   &
            dhdt,rho,t,po,us,vs,tn,dm2d 
        REAL,    DIMENSION (its:ite,kts:kte,1:ens4),INTENT (INOUT)    ::   &
            omeg
        REAL,    DIMENSION (its:ite,kts:kte)       ,INTENT (INOUT)    ::   &
            q,qo
        REAL,    DIMENSION (its:ite)               ,INTENT (IN   )    ::   &
            ccn,Z1,PSUR,xland,xlons,xlats, h_sfc_flux,le_sfc_flux,tsur,dx,  &
            stochastic_sig 
        REAL,    DIMENSION (its:ite)               ,INTENT (INOUT)    ::   &
            zws,ztexec,zqexec
        REAL                                       ,INTENT (IN   )    ::   &
            dtime,entr_rate_input 
        !
        ! local ensemble dependent variables in this routine
        real,    dimension (its:ite,1:maxens2) ::                         &
            edtc
        real,    dimension (its:ite,1:ensdim)        ::           &
            xf_ens,pr_ens
        !
        !*******the following are your basic environmental
        !          variables. They carry a "_cup" if they are
        !          on model cloud levels (staggered). They carry
        !          an "o"-ending (z becomes zo), if they are the forced
        !          variables. They are preceded by x (z becomes xz)
        !          to indicate modification by some typ of cloud
        !
        ! z                 = heights of model levels
        ! q           = environmental mixing ratio
        ! qes         = environmental saturation mixing ratio
        ! t           = environmental temp
        ! p           = environmental pressure
        ! he          = environmental moist static energy
        ! hes         = environmental saturation moist static energy
        ! z_cup       = heights of model cloud levels
        ! q_cup       = environmental q on model cloud levels
        ! qes_cup     = saturation q on model cloud levels
        ! t_cup       = temperature (Kelvin) on model cloud levels
        ! p_cup       = environmental pressure
        ! he_cup = moist static energy on model cloud levels
        ! hes_cup = saturation moist static energy on model cloud levels
        ! gamma_cup = gamma on model cloud levels
        !
        !
        ! hcd = moist static energy in downdraft
        ! zd normalized downdraft mass flux
        ! dby = buoancy term
        ! entr = entrainment rate
        ! zd   = downdraft normalized mass flux
        ! entr= entrainment rate
        ! hcd = h in model cloud
        ! bu = buoancy term
        ! zd = normalized downdraft mass flux
        ! gamma_cup = gamma on model cloud levels
        ! qcd = cloud q (including liquid water) after entrainment
        ! qrch = saturation q in cloud
        ! pwd = evaporate at that level
        ! pwev = total normalized integrated evaoprate (I2)
        ! entr= entrainment rate
        ! z1 = terrain elevation
        ! entr = downdraft entrainment rate
        ! jmin = downdraft originating level
        ! kdet = level above ground where downdraft start detraining
        ! psur        = surface pressure
        ! z1          = terrain elevation
        ! pr_ens = precipitation ensemble
        ! xf_ens = mass flux ensembles
        ! massfln = downdraft mass flux ensembles used in next timestep
        ! omeg = omega from large scale model
        ! mconv = moisture convergence from large scale model
        ! zd      = downdraft normalized mass flux
        ! zu      = updraft normalized mass flux
        ! dir     = "storm motion"
        ! mbdt    = arbitrary numerical parameter
        ! dtime   = dt over which forcing is applied
        ! iact_gr_old = flag to tell where convection was active
        ! kbcon       = LFC of parcel from k22
        ! k22         = updraft originating level
        ! ichoice       = flag if only want one closure (usually set to zero!)
        ! dby = buoancy term
        ! ktop = cloud top (output)
        ! xmb    = total base mass flux
        ! hc = cloud moist static energy
        ! hkb = moist static energy at originating level
        logical :: keep_going
        character*128 :: ierrc(its:ite)

        real,    dimension (its:ite,kts:kte) ::                           &
            entr_rate_2d,mentrd_rate_2d                                    &
            , he, hes, qes, z,heo,heso,qeso,zo, zu,zd                       &
            ,xhe,xhes,xqes,xz,xt,xq                                         &
            , qes_cup, q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup &
            ,qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,       gammao_cup,tn_cup &
            ,xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup                        &
            ,xt_cup,hcot,                                                   &

            dby,hc,clw_all,                                                &
            dbyo,qco,qrcdo,hcdo,qcdo,dbydo,hco,                            &
            xdby,xzu,xzd,   xhc, cupclw,pwo_eff,                           &

            ! cd  = detrainment function for updraft
            ! cdd = detrainment function for downdraft
            ! dellat = change of temperature per unit mass flux of cloud ensemble
            ! dellaq = change of q per unit mass flux of cloud ensemble
            ! dellaqc = change of qc per unit mass flux of cloud ensemble

            cd,cdd,dellah,dellaq,dellat,dellaqc,dsubq,dsubh,                &
            u_cup,v_cup,uc,vc,ucd,vcd,dellu,dellv,                          &
            up_massentr,up_massdetr,dd_massentr,dd_massdetr

            ! aa0 cloud work function for downdraft
            ! edt = epsilon
            ! aa0     = cloud work function without forcing effects
            ! aa1     = cloud work function with forcing effects
            ! xaa0    = cloud work function with cloud effects
            ! edt     = epsilon

        real,    dimension (its:ite) ::                            &
            edt,edtx,aa1,aa0,xaa0,hkb,                               &
            hkbo,xhkb,qkb, pwevo,bu,bud,cap_max,xland1,              &
            cap_max_increment,psum,psumh,sigd,mconv,rescale_entrain

        integer,    dimension (its:ite) ::                         &
            kzdown,kdet,kb, ierr2,ierr3,kbmax
        integer :: iloop,nall,iedt,nens,nens3,ki,I,K,KK,iresult,nvar,nvarbegin
        integer :: jprnt,k1,k2,kbegzu,kdefi,kfinalzu,kstart,jmini,imid,k_free_trop

        real :: day,dz,dzo,radius,entrd_rate,mentrd_rate,             &
            zcutdown,depth_min,zkbmax,z_detr,zktop,               &
            massfld,dh,cap_maxs,trash,frh,xlamdd,radiusd,frhd,effec_entrain
        real :: detdo1,detdo2,entdo,dp,subin,detdo,entup,             &
            detup,subdown,entdoj,entupk,detupk,totmas
        real :: tot_time_hr,beta,entr_rate,wmeanx
        real :: dts,fpi,denom,denomU

        real,    dimension (its:ite,1:maxens3) ::  xff_mid
        real,    dimension (kts:kte)   :: dummy1,dummy2
        integer :: iversion,bl=1,fa=2,step
        real :: umean,T_star

        real, dimension (its:ite)         :: aa0_bl,aa1_bl,tau_bl,tau_ecmwf,wmean,aa1_fa,aa1_tmp,hkbo_x &
                                                ,aa2,aa3,cin0,cin1,edtmax,edtmin
        real, dimension (its:ite,kts:kte),intent(IN) :: tn_bl, qo_bl
        real, dimension (its:ite,kts:kte) :: tn_x, qo_x, qeso_x, heo_x, heso_x &
                                    ,qeso_cup_x,qo_cup_x, heo_cup_x,heso_cup_x&
                                    ,gammao_cup_x,tn_cup_x,hco_x,DBYo_x

        real, dimension (its:ite,kts:kte) :: xhe_x,xhes_x,xt_x,xq_x,xqes_x, &
                                    xqes_cup_x,xq_cup_x,xhe_cup_x,xhes_cup_x,gamma_cup_x,xt_cup_x
        real, dimension (its:ite)         ::xaa0_x,xk_x

        real, dimension(its:ite) :: xf_dicycle,mbdt
        real :: C_up, E_dn,G_rain,trash2, pgc,bl2dp
        character(len=2) :: cty
        !
        real   :: dsubh_aver,dellah_aver,x_add,cap_max_inc,tlll,plll,rlll,tlcl,plcl,dzlcl,zlll
        integer:: start_k22,start_level(its:ite)
        real,    dimension (its:ite,kts:kte) ::  dtempdz
        integer, dimension (its:ite,kts:kte) ::  k_inv_layers
        integer :: ipr=0,jpr=0,fase,shift

        real,    dimension (its:ite,kts:kte) ::  vvel2d,tempco,tempcdo
        real,    dimension (its:ite        ) ::  vvel1d

        real,    dimension (its:ite,kts:kte) ::  p_liq_ice,melting_layer,melting
        real,    dimension (its:ite,kts:kte) ::  c1d
        real,    dimension (its:ite,kts:kte) :: prec_flx,evap_flx

        real,    dimension (its:ite) :: lambau
        real,    dimension (its:ite,kts:kte) ::  up_massentru,up_massdetru,&
                                        dd_massentru,dd_massdetru

        !- atmos composition treatment
        real, dimension (mtp,its:ite,kts:kte),intent (in)    ::   se_chem
        real, dimension (mtp,its:ite,kts:kte),intent (inout) ::   out_chem

        !-locals
        integer :: ispc
        real, dimension (mtp,its:ite,kts:kte) ::   se_cup_chem,sc_up_chem,sc_dn_chem,pw_up_chem &
                                        ,pw_dn_chem 
        real, dimension (mtp) ::   min_tend_chem
        real, dimension (mtp,its:ite)         ::  tot_pw_up_chem,tot_pw_dn_chem 
        real, dimension (mtp,kts:kte)         ::  trcflx_in,sub_tend,ddtr,ddtr_upd,zenv_diff,fp_mtp,fm_mtp
        real, dimension (kts:kte)             ::  aa,bb,cc,ddu,ddv,ddh,ddq,fp,fm
        real, dimension (its:ite,kts:kte)     ::  massflx,zenv
        real                                  ::  evap_(mtp),wetdep_(mtp),trash_(mtp),trash2_(mtp) &
                                                    ,massi,massf,dtime_max,residu_(mtp),evap,wetdep
        !----------------------------------------------------------------------
        integer, dimension (its:ite)      ,intent (inout)  :: &
                ierr              &
                ,jmin              &
                ,klcl              &
                ,k22               &
                ,kbcon             &
                ,ktop              &
                ,kstabi            &
                ,kstabm

        real,  dimension (its:ite)        ,intent (inout)  :: &
                xmb               &
                ,edto              &
                ,pwavo
        real,  dimension (its:ite,kts:kte),intent (inout)  :: &
                po_cup            &
                ,up_massentro      &
                ,up_massdetro      &
                ,dd_massentro      &
                ,dd_massdetro      &
                ,zuo               &
                ,zdo               &
                ,pwo               &
                ,pwdo              &
                ,qrco              &
                ,tup               &
                ,clfrac
        !-- debug/diag
        real,  dimension (its:ite)        ,intent (inout)  :: &
                aa0_,aa1_,aa2_,aa3_,aa1_bl_,aa1_cin_,tau_bl_,tau_ec_
        real,  dimension (its:ite,kts:kte) :: dtdt,dqdt
        real :: s1,s2,q1,q2,rzenv
        real :: alp0,beta1,beta2,dp_p,dp_m,delt1,delt2
        integer :: status,istep,lstep
        !----------------------------------------------------------------------

        !--- maximum depth (mb) of capping inversion (larger cap = no convection)
        if(trim(cumulus) == 'deep') then
            cap_maxs=150.
            cap_max_inc=20.
        endif
        if(trim(cumulus) == 'mid' ) then
            cap_maxs=150.
            cap_max_inc= 10.
        endif
        do i=its,itf
            kbmax  (i) = 1
            kstabm (i) = ktf-1
            ierr2  (i) = 0
            ierr3  (i) = 0
            xland1 (i) = xland(i) ! 1.
            cap_max(i) = cap_maxs
            lambau (i) = 2.
            ierrc  (i) = "ierrtxt"
            aa0    (i) = 0.0
            aa1    (i) = 0.0
            aa2    (i) = 0.0
            aa3    (i) = 0.0
            aa1_bl (i) = 0.0
            aa1_fa (i) = 0.0
            aa0_bl (i) = 0.0
            cin1   (i) = 0.0
            xk_x   (i) = 0.0
            edt    (i) = 0.0
            edto   (i) = 0.0
            tau_bl (i) = 0.0
            tau_ecmwf (i) = 0.0
            xf_dicycle(i) = 0.0
            z     (i,:) = zo(i,:)
            xz    (i,:) = zo(i,:)
            hcdo  (i,:) = 0.0
            cupclw(i,:) = 0.0
            qrcdo (i,:) = 0.0
            hcot  (i,:) = 0.0
            c1d   (i,:) = 0.0
            xf_ens(i,:) = 0.0
            pr_ens(i,:) = 0.0
            cap_max_increment(i)=cap_max_inc
        enddo
        !
        !--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
        !    base mass flux
        !
        !-- note : to make the evaporation stronger => increase "edtmin"
        if(trim(cumulus) == 'mid' ) then
            !edtmin=0.2 ;  edtmax=0.4 ! P9_P6
            !edtmin=0.75;  edtmax=0.9 ! P10
            edtmin(:)=0.1 ;  edtmax(:)=0.9 ! orig
            if(c0_mid < 1.e-8) edtmin(:)=0.0
        endif
        if(trim(cumulus) == 'deep') then
            DO i=its,itf
                if(xland(i) > 0.99 ) then !- over water
                    !edtmin(i)=0.1 ;  edtmax(i)=0.2 ! P9_P6
                    !edtmin(i)=0.05;  edtmax(i)=0.2 ! P10
                    edtmin(i)=0.10;  edtmax(i)=0.9 ! orig
                else!- over land
                    !edtmin(i)=0.1 ;  edtmax(i)=0.3 ! P9_P6
                    !edtmin(i)=0.05;  edtmax(i)=0.2 ! P10
                    edtmin(i)=0.10;  edtmax(i)=0.9 ! orig
                endif
            ENDDO
        endif
        !
        !--- minimum depth (m), clouds must have
        !
        if(trim(cumulus) == 'deep') depth_min=1000.
        if(trim(cumulus) == 'mid' ) depth_min=500.
        !
        !--- max height(m) above ground where updraft air can originate
        !
        if(trim(cumulus) == 'deep') zkbmax=4000.
        if(trim(cumulus) == 'mid' ) zkbmax=3000.
        !
        !--- height(m) above which no downdrafts are allowed to originate
        !
        zcutdown=3000.
        !
        !--- depth(m) over which downdraft detrains all its mass
        !
        z_detr=1000. !P9_P6
        if(trim(cumulus) == 'deep') z_detr= 1000. !500. DNDP
        if(trim(cumulus) == 'mid' ) z_detr= 300.
        !
        !- mbdt ~ xmb * timescale
        do i=its,itf
            mbdt(i)= 0.1!*dtime*xmb_nm1(i)
            !mbdt(i)= 100.*(p_cup(i,kbcon(i))-p_cup(i,kbcon(i)+1))/(g*dtime)
            !mbdt(i)= 0.1*mbdt(i)
        enddo

        !
        !--- environmental conditions, FIRST HEIGHTS
        !--- calculate moist static energy, heights, qes
        !
        call cup_env(z ,qes ,he ,hes ,t ,q ,po,z1 ,psur,ierr,-1,itf,ktf,its,ite, kts,kte)
        call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, psur,ierr,-1,itf,ktf,its,ite, kts,kte)
        !
        !--- environmental values on cloud levels
        !
        call cup_env_clev(t,qes,q,he,hes,z,po,qes_cup,q_cup,he_cup, &
            us,vs,u_cup,v_cup,                                     &
            hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,tsur,         &
            ierr,z1,itf,ktf,its,ite, kts,kte, host_model)

        call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, heo_cup,   &
            us,vs,u_cup,v_cup,                                           &
            heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,tsur,          &
            ierr,z1,itf,ktf,its,ite, kts,kte, host_model)


        !IF( MAPL_AM_I_ROOT()) then
        !     print*,"1GF =============================================================="
        !     do k=kts,ktf-1
        !         write(15,10) k,po(1,k),po_cup(1,k),100.*(po_cup(1,k)-po_cup(1,k+1))/g,dm2d (1,k)
        !	       10 format(1x,I4,4E15.7)
        !     end do
        !     print*,"2GF =============================================================="
        !     flush(6)
        !ENDIF
        !
        !--- partition between liq/ice cloud contents
        call get_partition_liq_ice(ierr,tn,z1,zo_cup,po_cup,p_liq_ice,melting_layer,&
                            itf,ktf,its,ite,kts,kte,cumulus)
        !     
        do i=its,itf
            if(ierr(i).eq.0)then
                do k=kts,ktf
                    if(zo_cup(i,k).gt.zkbmax+z1(i))then
                        kbmax(i)=k ; go to 25
                    endif
                enddo
                25       continue
                !
                !--- level where detrainment for downdraft starts
                !
                do k=kts,ktf
                    if(zo_cup(i,k).gt.z_detr+z1(i))then
                        kdet(i)=k; go to 26
                    endif
                enddo
                26       continue
            endif
        enddo
        !
        !--- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
        !
        start_k22 = 2
        k22(:)=0
        DO 36 i=its,itf
            IF(ierr(i).eq.0)THEN
                k22(i)=maxloc(HEO_CUP(i,start_k22:kbmax(i)+1),1)+start_k22-1
                IF(k22(i) > kbmax(i))then
                    !- let's try k22=start_k22 for the cases k22>kbmax
                    k22(i)= start_k22
                    cycle
                endif
            endif
        36   CONTINUE
        !
        !------- DETERMINE LCL for the air parcels around K22
        !
        DO i=its,itf
            klcl(i) = k22(i) ! default value
            if(ierr(i) == 0)then
                !tlll, rlll,plll - temp, water vapor and pressure of the source air parcel
                x_add = float(use_excess)*zqexec(i)
                call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),q_cup (i,kts:kte),rlll,k22(i),x_add)
                x_add = float(use_excess)*ztexec(i)
                call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),t_cup (i,kts:kte),tlll,k22(i),x_add)
                call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),p_cup (i,kts:kte),plll,k22(i))

                call get_lcl(tlll,100.*plll,rlll,tlcl,plcl,dzlcl)
                !print*,"MID",tlll,100.*plll,rlll,tlcl,plcl,dzlcl; flush(6)

                !-get LCL index
                if(dzlcl >= 0.) then ! LCL found (if dzlcl<0 => not found)
                    call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),z_cup (i,kts:kte),zlll,k22(i))
                    loop0:       do k=kts,ktf
                        if(z_cup(i,k).gt.zlll+dzlcl)then
                            klcl(i)=max(k,k22(i))
                            exit loop0
                        endif
                    enddo loop0
                    klcl(i)=min(klcl(i),ktf-4)
                endif
            endif
            !write(12,111)'MDlcl',tlcl,plcl,dzlcl,klcl(i),ierr(i)
            111 format(1x,A5,3F10.2,2i4)
        Enddo
        !-- check if LCL is below PBL height for shallow convection
        if(USE_LCL .and. trim(cumulus) == 'mid' )then
            do i=its,itf
                if(ierr(i).eq.0)then
                    if(klcl(i) > max(1,kpbl(i)+1)) then
                        ierr(i)=21
                        ierrc(i)='for mid convection:  LCL height > PBL height'
                    endif
                endif
            ENDDO
        endif
        !
        !--- DETERMINE THE moist static energy of air parcels at source level
        !
        do i=its,itf
            if(ierr(I) /= 0)cycle
            x_add = float(use_excess)*(xlv*zqexec(i)+cp*ztexec(i))
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),he_cup (i,kts:kte),hkb (i),k22(i),x_add)
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup(i,kts:kte),hkbo(i),k22(i),x_add)
        enddo
        !
        !--- define entrainment/detrainment profiles for updrafts
        !
        !- initial entrainment/detrainment
        entr_rate         = entr_rate_input
        entr_rate_2d(:,:) = entr_rate_input
        cd          (:,:) = entr_rate_input
        l_SIG:DO fase = 1,2
            IF(fase == 1) THEN
                !- default value for entrainmente re-scale
                rescale_entrain(:)=1.
                DO i=its,itf
                    IF(ierr(i) /= 0) CYCLE
                    do k=kts,ktf

                        frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
                        !-------------------------------------------
                        if(ENTRNEW) then
                            !- v 2
                            if(k >= klcl(i)) then
                                    entr_rate_2d(i,k)=entr_rate*(1.3-frh)*(qeso_cup(i,k)/qeso_cup(i,klcl(i)))**3
                            else
                                    entr_rate_2d(i,k)=entr_rate*(1.3-frh)
                            endif
                            cd(i,k)=0.75e-4*(1.6-frh)
                        else
                            !- v 1
                            entr_rate_2d(i,k)=max(entr_rate*(1.3-frh)*max(min(1.,(qeso_cup(i,k)/qeso_cup(i,klcl(i)))**1.25),0.1),1.e-5)
                            if(trim(cumulus) == 'deep') cd(i,k)=1.e-2*entr_rate
                            if(trim(cumulus) == 'mid' ) cd(i,k)=0.75*entr_rate_2d(i,k)
                        endif
                    enddo
                ENDDO
            ELSE  !- for 2nd fase (re-scale entrainemnt for the scale dependence approach, if required)
                DO i=its,itf
                    if(ierr(i) /= 0 .or. rescale_entrain(i) < 1.0001) cycle
                    !-rescale entr/detrain associated with the scale dependence approach (see below)
                    entr_rate_2d(i,kts:ktf) = rescale_entrain(i)*entr_rate_2d(i,kts:ktf)
                    cd          (i,kts:ktf) = rescale_entrain(i)*cd          (i,kts:ktf)
                enddo
            ENDIF
            !
            !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
            !
            CALL cup_cloud_limits(cumulus,ierrc,ierr,cap_max_increment,cap_max,entr_rate  &
                            ,heo_cup,heso_cup,qo_cup,qeso_cup,po,po_cup,zo_cup,heo,hkbo &
                            ,entr_rate_2d,hcot,k22,kbmax,klcl,kbcon,ktop             &
                            ,use_excess,zqexec,ztexec,xland,itf,ktf,its,ite, kts,kte)

            if( USE_SCALE_DEP == 0) then
                sig(:)=1.; EXIT l_SIG
            endif

            !
            !--- SCALE DEPENDENCE FACTOR (SIG), version new
            !
            if( USE_SCALE_DEP == 1 ) then
                do i=its,itf
                    sig(i) = 0.
                    if(ierr(i) /= 0) cycle
                    sig(i)= (1.0-0.9839*exp(-0.09835*(dx(i)/1000.)))
                    if (stochastic_sig(i) /= 1.0) then
                        sig(i) = sig(i)**(stochastic_sig(i)*MAX(1.0,2.5*sig(i)))
                    endif
                    sig(i)= max(0.001,min(sig(i),1.))
                    !print*,"FORM2=",sig(i),dx(i)
                enddo
                EXIT l_SIG
            endif

            if(fase == 2) EXIT l_SIG
            !
            !--- SCALE DEPENDENCE FACTOR (SIG), version original
            !
            if( USE_SCALE_DEP == 2 ) then
                !- get the effective entraiment between klcl and kbcon
                do i=its,itf
                    sig(i) = 1.
                    IF(ierr(i) /= 0)cycle
                    effec_entrain=0.0
                    do k=klcl(i),kbcon(i)
                        dp = po_cup(i,k)-po_cup(i,k+1)
                        effec_entrain=effec_entrain+entr_rate_2d(i,k)*dp
                    enddo
                    effec_entrain=effec_entrain/(po_cup(i,klcl(i))-po_cup(i,kbcon(i)+1))

                    !- scale dependence factor
                    radius =0.2/effec_entrain
                    frh    =3.14*radius*radius/(dx(i)*dx(i))
                    !  print*,"frh1=",frh,effec_entrain,(1.-frh)**2
                    if(frh .gt. 0.7)then
                        frh=0.7
                        radius=sqrt(frh*dx(i)*dx(i)/3.14)! dx*sqrt(frh/3.14)
                        rescale_entrain(i) =(0.2/radius)/effec_entrain

                        !print*,"frh2=",frh,effec_entrain,rescale_entrain(i)
                    endif
                    !- scale dependence factor
                    sig (i)=(1.-frh)**2 
                    !print*,"FORM1=",sig(i),dx(i)
                enddo
            endif
        ENDDO l_SIG ! fase loop


        !--- define entrainment/detrainment profiles for downdrafts
        if(ENTRNEW) then
            mentrd_rate = entr_rate*0.3
        else
            mentrd_rate = entr_rate
        endif
        cdd(:,:) = mentrd_rate
        !- scale dependence factor
        sigd(:)  = 1.
        if( .not. DOWNDRAFT) sigd(:) = 0.0


        !--- update hkb/hkbo in case of k22 is redefined in 'cup_kbon'
        do i=its,itf
            IF(ierr(i) /= 0) CYCLE
            x_add = float(use_excess)*(xlv*zqexec(i)+cp*ztexec(i))
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),he_cup (i,kts:kte),hkb (i),k22(i),x_add)
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup(i,kts:kte),hkbo(i),k22(i),x_add)
        enddo
        !
        !--- increase detrainment in stable layers
        !
        CALL cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,itf,ktf,its,ite, kts,kte)
        !
        !--- get KTOP for mid and deep plumes
        !
        CALL rates_up_pdf(cumulus,ktop,ierr,po_cup,entr_rate_2d,hkbo,heo,heso_cup,zo_cup, &
                kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl,klcl,hcot)
        !
        !--- option for using the inversion layers as a barrier for the convection development
        !--- (only for mid congestus)
        !
        IF(trim(cumulus) == 'mid') THEN
            if(USE_INV_LAYERS) then
                !
                !--- get inversion layers
                call get_inversion_layers(cumulus,ierr,psur,po_cup,tn_cup,zo_cup,k_inv_layers,&
                                    dtempdz,itf,ktf,its,ite, kts,kte)
                DO i=its,itf
                    if(ierr(i) /= 0)cycle
                    ktop(i) = min(ktop(i),k_inv_layers(i,mid))
                    !print*,"ktop=",ktop(i),k_inv_layers(i,mid);flush(6)
                enddo
            endif
            !
            !-- check if ktop is above 450hPa layer for mid convection
            !
            do i=its,itf
                if(ierr(i) /= 0)cycle
                !print*,"sta=",Kbcon(i),kstabm(i),kstabi(i),p_cup(i,ktop(i)),z_cup(i,kstabi(i))
                if(po_cup(i,ktop(i)) < 450.) then
                    ierr(i)=25
                    ierrc(i)='mid convection with cloud top above 450 hPa'
                endif
            ENDDO

        ENDIF
        !
        !-- last checks for ktop (deep and mid)
        !
        DO i=its,itf
            if(ierr(i) /= 0)cycle
            if(ktop(i) < kbcon(i)+5)then
                ierr(i)=5
                ierrc(i)='ktop too small'
            endif
        ENDDO

        DO i=its,itf
            if(ierr(i) /= 0)cycle
            if(po_cup(i,ktop(i)) > 750.) then
                ierr(i)=55
                ierrc(i)='ktop too low for deep/mid'
            endif
        ENDDO
        !
        !-- avoid double-counting with shallow scheme (deep and mid)
        !
        DO i=its,itf
            if(ierr(i) /= 0)cycle
            if(last_ierr(i) == 0) then
                !--- if 'mid' => last was 'shallow'
                ! if(trim(cumulus) == 'mid' .and. po_cup(i,ktop(i)) > 700.) then
                !   ierr(i)=27
                !   ierrc(i)='avoiding double-counting shallow and mid'
                ! endif
                !--- if 'mid' => last was 'shallow'
                if(trim(cumulus) == 'mid' ) then
                    ierr(i)=27
                    ierrc(i)='avoiding double-counting deep and mid'
                endif
            endif
        ENDDO
        !
        !--- determine the normalized mass flux profile for updraft
        !
        do i=its,itf
            zuo(i,:)=0.
            if(ierr(i) /= 0) cycle
            ! CALL get_zu_zd_pdf(trim(cumulus),trim(cumulus)//"_up",ierr(i),k22(i),ktop(i),zuo(i,kts:kte),kts,kte,ktf  &
            !             ,kpbl(i),k22(i),kbcon(i),klcl(i),po_cup(i,kts:kte),psur(i),xland(i))
        enddo

        do i=its,itf
            if(ierr(i) /= 0) cycle
            xzu(i,:)= zuo(i,:)
            zu (i,:)= zuo(i,:)
        enddo
        !
        ! calculate mass entrainment and detrainment
        !
        CALL get_lateral_massflux(itf,ktf, its,ite, kts,kte                   &
                        ,ierr,ktop,zo_cup,zuo,cd,entr_rate_2d        &
                        ,up_massentro, up_massdetro ,up_massentr, up_massdetr &
                        ,cumulus,kbcon,k22,up_massentru,up_massdetru,lambau)
        !--- start_level
        if(trim(cumulus) == 'deep') start_level(:)=  KLCL(:)!k22(:)
        if(trim(cumulus) == 'mid' ) start_level(:)=  KLCL(:)!kbcon(:)

        uc  =0.
        vc  =0.
        hc  =0.
        hco =0.
        do i=its,itf
            IF(ierr(i).eq.0)THEN
                do k=kts,start_level(i)
                    hc (i,k) =hkb (i)
                    hco(i,k) =hkbo(i)
                    !-get uc and vc as average between layers below k22
                    call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),u_cup(i,kts:kte),uc(i,k),k22(i))
                    call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),v_cup(i,kts:kte),vc(i,k),k22(i))
                enddo
            ENDIF
        enddo
        !
        !--- 1st guess for moist static energy and dbyo (not including ice phase)
        !
        do i=its,itf
            if(ierr(i) /= 0) cycle
            do k=start_level(i)  +1,ktop(i) + 1  ! mass cons option
                denom=(zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
                if(denom > 0.0)then
                    hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                                                    up_massentro(i,k-1)*heo(i,k-1))  / denom
                else
                    hco(i,k)= hco(i,k-1)
                endif
            enddo
            do k=ktop(i)+2,ktf
                hco (i,k)=heso_cup(i,k)!=heo_cup(i,k)
            enddo
        enddo
        !
        !--- Get bouyancy of updrafts
        !
        call get_bouyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
            ,hco,heo_cup,heso_cup,dbyo,zo_cup)

        !--- get "c1d" profile ----------------------------------------
        if(trim(cumulus) == 'deep' .and. USE_C1D) then
            Do i=its,itf
                if(ierr(i) /= 0) cycle
                c1d(i,kbcon(i)+1:ktop(i)-1)=c1
            ENDDO
        endif
        !----------------------------------------------------------------

        !--- calculate moisture properties of updraft
        IF(CLOUDMP==1) then
            call cup_up_moisture(cumulus,klcl,ierr,zo_cup,qco,qrco,pwo,pwavo,hco,tempco,xland,  &
                po,p_cup,kbcon,ktop,cd,dbyo,clw_all,                             &
                t_cup,qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,                     &
                ZQEXEC,use_excess,ccn,rho,up_massentr,up_massdetr,psum,psumh,c1d,&
                1,itf,ktf,ipr,jpr,its,ite, kts,kte)

        ELSEIF(CLOUDMP==2) then
            !.. do i=its,itf
            !.. qco(i,:) = qo_cup(i,:)
            !.. tempco(i,:) = t_cup(i,:)
            !.. vvel2d(i,:) = 0.
            !.. if(ierr(i) /= 0) cycle
            !.. start_level(i)=klcl(i)
            !.. x_add =float(use_excess)*zqexec(i)
            !.. call get_cloud_bc(kts,kte,qo_cup (i,kts:kte),qco(i,kts),k22(i),x_add)

            !.. do k=kts,start_level(i)
            !.. qco(i,k) = qco(i,kts) !- air parcel water vapor
            !.. tempco(i,k) = (1./cp)*(hco(i,k)-g*zo_cup(i,k)-xlv*qco(i,k))
            !.. vvel2d(i,k) = zws(i)
            !.. !---- temporary !<<<<<<<<<<<<<<<<<<<<<
            !.. !tempco(i,k)=tempco(i,k)+2.
            !.. print*,"Ks=",ktop(i),kbcon(i),start_level(i)
            !.. print*,"CL=",k, qco(i,k),tempco(i,k),hco(i,k),zo_cup(i,k)
            !.. !---- temporary !<<<<<<<<<<<<<<<<<<<<<
            !.. enddo

            !.. ENDDO
            !.. call GTMP_2_GFCONVPAR_interface(its,ite,itf, kts,kte,ktf,cumulus &
            !.. ,ierr,klcl,kbcon,ktop,k22,kpbl             &
            !.. ,zo_cup,po_cup,tn_cup,qo_cup                 &
            !.. ,zo,tempco,qo,qco,qrco,hco,pwo,pwavo,rho       &
            !.. ,ztexec,zqexec,use_excess,zuo,up_massentr,up_massdetr,entr_rate_2d,cd &
                    !.. ,dbyo,GAMMAo_cup,qeso_cup,vvel2d,vvel1d)
            !.. RETURN
            !.. stop 333

        ENDIF
        do i=its,itf
            if(ierr(i).eq.0)then
                do k=kts,ktop(i)+1
                    cupclw(i,k)=qrco(i,k)
                enddo
            endif
        enddo
        !
        !--- get melting profile
        call get_melting_profile(ierr,tn_cup,po_cup, p_liq_ice,melting_layer,qrco    &
                        ,pwo,edto,pwdo,melting                               &
                        ,itf,ktf,its,ite, kts,kte, cumulus                   )
        !---meltglac-----------------------------------------------------
        !
        !
        do i=its,itf
            if(ierr(i) /= 0) cycle
            do k=start_level(i)  +1,ktop(i)+1  ! mass cons option
                denom =(zu(i,k-1)-.5*up_massdetr (i,k-1)+up_massentr (i,k-1))
                denomU=(zu(i,k-1)-.5*up_massdetru(i,k-1)+up_massentru(i,k-1))
                if(denom > 0.0 .and. denomU >0.0)then

                    hc (i,k)=(hc (i,k-1)*zu (i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1) + &
                                                    up_massentr(i,k-1)*he(i,k-1))/ denom

                    hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                                                    up_massentro(i,k-1)*heo(i,k-1))  / denom
                                !assuming zuo=zu,up_massdetro=up_massdetr, ...
                                !(zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))

                    uc(i,k)=(uc(i,k-1)*zu(i,k-1)-.5*up_massdetru(i,k-1)*uc(i,k-1)+     &
                                                up_massentru(i,k-1)*us(i,k-1)      &
                        -pgcon*.5*(zu(i,k)+zu(i,k-1))*(u_cup(i,k)-u_cup(i,k-1))) /  denomU

                    vc(i,k)=(vc(i,k-1)*zu(i,k-1)-.5*up_massdetru(i,k-1)*vc(i,k-1)+     &
                                                up_massentru(i,k-1)*vs(i,k-1)      &
                        -pgcon*.5*(zu(i,k)+zu(i,k-1))*(v_cup(i,k)-v_cup(i,k-1))) /  denomU


                else
                    hc (i,k)= hc (i,k-1)
                    hco(i,k)= hco(i,k-1)
                    uc (i,k)= uc (i,k-1)
                    vc (i,k)= vc (i,k-1)
                endif
                !---meltglac-------------------------------------------------
                !- includes glaciation effects on HC,HCO
                !                    ------ ice content --------
                !print*,"H=",hc (i,k),(1.-p_liq_ice(i,k))*qrco(i,k)*xlf,hc (i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf
                hc (i,k)= hc (i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf
                hco(i,k)= hco(i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf

            enddo
            !
            do k=ktop(i)+2,ktf
                hc  (i,k)= hes_cup(i,k)!= he_cup(i,k)
                uc  (i,k)=   u_cup(i,k)
                vc  (i,k)=   v_cup(i,k)
                hco (i,k)=heso_cup(i,k)!=heo_cup(i,k)
                zu  (i,k)=0.
                zuo (i,k)=0.
            enddo
        enddo
        !
        !--- Get bouyancy of updrafts
        !
        call get_bouyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
            ,hc,he_cup,hes_cup,dby,z_cup)
        call get_bouyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
            ,hco,heo_cup,heso_cup,dbyo,zo_cup)
        !
        !
        !--- calculate in-cloud/updraft air temperature for vertical velocity
        !
        do i=its,itf
            if(ierr(i) == 0)then
                do k=kts,ktf
                    tempco (i,k) = (1./cp)*(hco (i,k)-g*zo_cup(i,k)-xlv*qco (i,k))
                enddo
            else
                tempco (i,:)=tn_cup(i,:)
            endif
        enddo
        !
        !--- vertical velocity
        !
        call cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd,zo,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup &
                ,tempco,qco,qrco,qo,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

        !
        !--- DOWNDRAFT section
        !
        DO 37 i=its,itf
            kzdown(i)=0
            if(ierr(i).eq.0)then
                zktop=(zo_cup(i,ktop(i))-z1(i))*.6
                zktop=min(zktop+z1(i),zcutdown+z1(i))
                do k=kts,ktf
                    if(zo_cup(i,k).gt.zktop)then
                        kzdown(i)=k
                        go to 37
                    endif
                enddo
            endif
        37   CONTINUE
        !
        !--- DOWNDRAFT ORIGINATING LEVEL - JMIN
        !
        call cup_minimi(heso_cup,k22,kzdown,jmin,ierr,itf,ktf,its,ite, kts,kte)

        DO 100 i=its,ite
            if(i.eq.ipr)then
                write(0,*)i,'k,p_cup(i,k),heo_cup(i,k),heso_cup(i,k),hco(i,k)'
                do k=kts,min(40,ktf)
                    write(0,*)k,p_cup(i,k),heo_cup(i,k),heso_cup(i,k),hco(i,k)
                Enddo
            endif
            IF(ierr(I).eq.0)THEN

                !----new jmin
                if(trim(cumulus) == 'deep' .and. MELT_GLAC) jmin(i)=max(jmin(i),maxloc(melting_layer(i,:),1))

                !----new jmin

                !
                !--- check whether it would have buoyancy, if there where
                !--- no entrainment/detrainment
                !
                jmini = jmin(i)
                keep_going = .TRUE.
                do while ( keep_going )
                    keep_going = .FALSE.
                    if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
                    if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
                    ki = jmini
                    hcdo(i,ki)=heso_cup(i,ki)
                    DZ=Zo_cup(i,Ki+1)-Zo_cup(i,Ki)
                    dh=0.
                    do k=ki-1,1,-1
                        hcdo(i,k)=heso_cup(i,jmini)
                        DZ=Zo_cup(i,K+1)-Zo_cup(i,K)
                        dh=dh+dz*(HCDo(i,K)-heso_cup(i,k))
                        if(dh.gt.0.)then
                            jmini=jmini-1
                            if ( jmini .gt. 5 ) then
                                keep_going = .TRUE.
                            else
                                ierr(i)=9
                                ierrc(i) = "could not find jmini9"
                                exit
                            endif
                        endif
                    enddo
                enddo
                jmin(i) = jmini
                if ( jmini .le. 5 ) then
                    ierr(i)=4
                    ierrc(i) = "could not find jmini4"
                endif
            ENDIF
        100   continue
        !
        ! - Must have at least depth_min m between cloud convective base
        !     and cloud top.
        !
        do i=its,itf
            IF(ierr(I).eq.0)THEN
                if ( jmin(i) - 1 .lt. kdet(i)) kdet(i) = jmin(i)-1
                if (-zo_cup(i,kbcon(i))+zo_cup(i,ktop(i)).lt.depth_min)then
                    ierr(i)=6
                    ierrc(i)="cloud depth very shallow"
                endif
            endif
        enddo

        if(trim(cumulus) == 'deep') beta=0.05
        if(trim(cumulus) == 'mid' ) beta=0.02
        !
        !--- this calls routine to get downdrafts normalized mass flux
        !
        do i=its,itf
            zd(i,:)=0.
            IF(ierr(i) /= 0) cycle
            ! call get_zu_zd_pdf(trim(cumulus),"DOWN",ierr(i),kdet(i),jmin(i),zdo(i,:),kts,kte,ktf&
            !             ,kpbl(i),k22(i),kbcon(i),klcl(i),po_cup(i,kts:kte),psur(i),xland(i))
        enddo
        !
        !--- this calls routine to get lateral mass fluxes associated with downdrafts
        !
        CALL get_lateral_massflux_down(itf,ktf, its,ite, kts,kte &
                                ,ierr,jmin,zo_cup,zdo,xzd,zd,cdd,mentrd_rate_2d      &
                                ,dd_massentro,dd_massdetro ,dd_massentr, dd_massdetr &
                                ,cumulus,mentrd_rate,dd_massentru,dd_massdetru,lambau)
        !
        !--- downdraft moist static energy + moisture budget
        !
        do i=its,itf
            hcdo (i,:)= heso_cup(i,:)
            ucd  (i,:)=    u_cup(i,:)
            vcd  (i,:)=    v_cup(i,:)
            dbydo(i,:)= 0.
        enddo

        do i=its,itf
            bud(i)=0.
            if(ierr(i)/= 0)cycle

            dbydo(i,jmin(i))=hcdo(i,jmin(i))-heso_cup(i,jmin(i))
            bud(i)=dbydo(i,jmin(i))*(zo_cup(i,jmin(i)+1)-zo_cup(i,jmin(i)))

            do ki=jmin(i)  ,kts,-1!do ki=jmin(i)-1,1,-1
                denom = zdo(i,ki+1)-0.5*dd_massdetro(i,ki)+dd_massentro(i,ki)
                denomU= zdo(i,ki+1)-0.5*dd_massdetru(i,ki)+dd_massentru(i,ki)
                !-tmp fix for denominator being zero
                if(denom > 0.0 .and. denomU >0.0)then
                    dzo=zo_cup(i,ki+1)-zo_cup(i,ki)

                    ucd(i,ki)=(ucd(i,ki+1)*zdo(i,ki+1)-.5*dd_massdetru(i,ki)*ucd(i,ki+1)+ &
                                                            dd_massentru(i,ki)*us (i,ki)    &
                                -pgcon*zdo(i,ki+1)*(us(i,ki+1)-us(i,ki)))   /  denomU
                    vcd(i,ki)=(vcd(i,ki+1)*zdo(i,ki+1)-.5*dd_massdetru(i,ki)*vcd(i,ki+1)+        &
                                                            dd_massentru(i,ki)*vs (i,ki)           &
                                -pgcon*zdo(i,ki+1)*(vs(i,ki+1)-vs(i,ki)))   /  denomU

                    hcdo(i,ki)=(hcdo(i,ki+1)*zdo(i,ki+1)-.5*dd_massdetro(i,ki)*hcdo(i,ki+1)+     &
                                                            dd_massentro(i,ki)*heo(i,ki))  /denom

                    dbydo(i,ki)=hcdo(i,ki)-heso_cup(i,ki)
                    !if(i.eq.ipr)write(0,*)'ki,bud = ',ki,bud(i),hcdo(i,ki)
                    bud(i)=bud(i)+dbydo(i,ki)*dzo
                else
                    ucd (i,ki)= ucd(i,ki+1)
                    vcd (i,ki)= vcd(i,ki+1)
                    hcdo(i,ki)=hcdo(i,ki+1)
                endif
            enddo
            if(bud(i).gt.0)then
                ierr(i)=7
                ierrc(i)='downdraft is not negatively buoyant '
            endif
        enddo
        !
        !--- calculate moisture properties of downdraft
        !
        call cup_dd_moisture(cumulus,ierrc,zdo,hcdo,heso_cup,qcdo,qeso_cup, &
            pwdo,qo_cup,zo_cup,dd_massentro,dd_massdetro,jmin,ierr,gammao_cup, &
            pwevo,pwavo,bu,qrcdo,qo,heo,tn_cup,1, itf,ktf,its,ite, kts,kte)
        !
        !
        !--- calculate workfunctions for updrafts
        !
        call cup_up_aa0(aa0,z_cup ,zu ,dby  ,GAMMA_CUP    ,t_cup   ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)
        call cup_up_aa0(aa1,zo_cup,zuo,dbyo ,GAMMAo_CUP   ,tn_cup  ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

        do i=its,itf
            if(ierr(i) /= 0) cycle
            if(aa1(i).eq.0.)then
                ierr(i)=17
                ierrc(i)="cloud work function zero"
            endif
        enddo
        !
        !--- calculate CIN for updrafts
        !
        call cup_up_aa0(cin0,z_cup ,zu ,dby  ,GAMMA_CUP    ,t_cup   ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,'CIN')
        call cup_up_aa0(cin1,zo_cup,zuo,dbyo ,GAMMAo_CUP   ,tn_cup  ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,'CIN')
        !
        !----trigger function for KE+CIN < 0 => no convection
        !
        IF( DICYCLE>1) THEN
            do i=its,itf
                if(ierr(i) /= 0) cycle
                print*,"cin=",cin0(i),0.5*zws(i)**2, omeg(i,kpbl(i),1)/(-g*rho(i,kpbl(i))) ;flush(6)
                if(cin0(i) + 0.5*zws(i)**2 < 0.)then !think about including the grid scale vertical velocity at KE calculation
                    ierr(i)=19
                    ierrc(i)="CIN negat"
                endif
            enddo
        ENDIF
        !
        !
        !--- calculate in-cloud/updraft and downdraft air temperature for vertical velocity
        !
        do i=its,itf
            if(ierr(i) == 0)then
                do k=kts,ktf
                    tempcdo(i,k) = (1./cp)*(hcdo(i,k)-g*zo_cup(i,k)-xlv*qcdo(i,k))
                enddo
            else
                tempcdo(i,:)=tn_cup(i,:)
            endif
        enddo
        !
        !--- diurnal cycle section
        !
        if(trim(cumulus)=='deep') then
            wmeanx=3.   ! m/s ! in the future change for Wmean == integral( W dz) / cloud_depth
            T_star=5.  ! T_star = temp scale in original paper = 1 K
        endif
        if(trim(cumulus)=='mid') then
            wmeanx=3!.4.
            T_star=40.!1
        endif

        !
        ! Bechtold et al 2008 time-scale of cape removal
        DO i=its,itf
            if(ierr(i) /= 0) cycle
            !- mean vertical velocity
            wmean(i) = wmeanx ! m/s ! in the future change for Wmean == integral( W dz) / cloud_depth
            !wmean(i) = min(max(vvel1d(i),1.),15.)
            !- time-scale cape removal from Bechtold et al. 2008
            tau_ecmwf(i)=( zo_cup(i,ktop(i))- zo_cup(i,kbcon(i)) ) / wmean(i)

            if(trim(cumulus)=='deep') tau_ecmwf(i)=max(tau_ecmwf(i),tau_deep)
            if(trim(cumulus)=='mid' ) tau_ecmwf(i)=max(tau_ecmwf(i),tau_mid)

            tau_ecmwf(i)= tau_ecmwf(i) * (1. + 1.66 * (dx(i)/(125*1000.)))! dx must be in meters
        enddo
        !
        DO i=its,itf
            if(ierr(i) /= 0) cycle
            if(xland(i) > 0.99 ) then !- over water
                umean= 2.0+sqrt(0.5*(US(i,1)**2+VS(i,1)**2+US(i,kbcon(i))**2+VS(i,kbcon(i))**2))
                tau_bl(i) = (zo_cup(i,kbcon(i))- z1(i)) /umean
            else !- over land
                tau_bl(i) =( zo_cup(i,ktop(i))- zo_cup(i,kbcon(i)) ) / wmean(i)
                !tau_bl(i)=max(tau_bl(i),1800.)
            endif
        ENDDO

        IF( (DICYCLE==1 .or. DICYCLE == 6) .or. (DICYCLE==0 .and. trim(cumulus)=='mid') ) THEN
            iversion=0
            !-- calculate pcape from BL forcing only
            call cup_up_aa1bl(iversion,aa1_bl,aa1_fa,aa1,t,tn,q,qo,dtime,po_cup,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup, &
                        rho,klcl,kpbl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,                          &
                        xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux,tau_bl,tau_ecmwf,t_star,cumulus,tn_bl,qo_bl  )
            if(DICYCLE==6) then
                call cup_up_aa1bl(iversion,aa0_bl,aa1_fa,aa1,t,tn_bl,q,qo_bl,dtime,po_cup,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup, &
                                rho,klcl,kpbl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,&
                                xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux,tau_bl,tau_ecmwf,t_star,cumulus,tn_bl,qo_bl   )
                where(ierr == 0) aa1_bl=0.50*aa1_bl+0.50*aa0_bl !mix1
            endif
            DO i=its,itf
                if(ierr(i) /= 0) cycle
                aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i)
                !           aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i) - cin1(i)
            ENDDO
            !==========================================================
        ELSEIF( DICYCLE==4) then!zhang(2002)

            !- T and Q profiles modified only by RAD+ADV tendencies
            DO i=its,itf
                if ( ierr(i) /= 0 ) CYCLE
                tn_x(i,kts:ktf) = tn(i,kts:ktf)-tn_bl(i,kts:ktf)+t(i,kts:ktf)
                qo_x(i,kts:ktf) = qo(i,kts:ktf)-qo_bl(i,kts:ktf)+q(i,kts:ktf)
            ENDDO

            !--- calculate moist static energy, heights, qes, ... only by free troposphere tendencies
            call cup_env(zo,qeso_x,heo_x,heso_x,tn_x,qo_x,po,z1, &
                psur,ierr,-1,itf,ktf, its,ite, kts,kte)
            !--- environmental values on cloud levels only by FT tendencies
            call cup_env_clev(tn_x,qeso_x,qo_x,heo_x,heso_x,zo,po,qeso_cup_x,qo_cup_x,heo_cup_x,           &
                        us,vs,u_cup,v_cup,                                    &
                        heso_cup_x,zo_cup,po_cup,gammao_cup_x,tn_cup_x,psur,tsur,  &
                        ierr,z1,itf,ktf,its,ite, kts,kte, host_model)
            !--- this is (DT_ve/Dt)_adv+rad
            DO i=its,itf
                IF(ierr(i) /= 0) cycle
                aa3(i)=0.
                DO k=max(kbcon(i),kts+1),ktop(i)
                    dp=- (log(100.*po(i,k))-log(100.*po(i,k-1))) !no units
                    aa3(i)=aa3(i)- (tn_cup_x(i,k)*(1.+0.608*qo_cup_x(i,k)) - &
                                    t_cup   (i,k)*(1.+0.608*q_cup   (i,k)) ) *dp/dtime !units = K/s
                !print*,"tve=",k,aa3(i),tn_cup_x(i,k)*(1.+0.608*qo_cup_x(i,k)),&
                !               t_cup (i,k)*(1.+0.608*q_cup   (i,k)),dp
                ENDDO
            ENDDO
            DO i=its,itf
                if(ierr(i)/= 0)CYCLE
                !- this is (DCAPE_env/Dt)_adv+rad
                !aa1_bl(i) = -aa3(i)
                !- Zhang threshold:  65 J/kg/hour => 65/(Rd *3600)= 63 10^-6 K/s
                aa1_bl(i) = aa3(i)-(63.e-6)!*1.5
                !print*,"dcape_env=",aa3(i),aa1_bl(i)
                if(xland(i) > 0.90 ) aa1_bl(i)=1.4*aa1_bl(i) !- over water
            ENDDO
            !--- this is (DT_ve/Dt)_cu
            DO i=its,itf
                dtdt(i,:)=0.
                dqdt(i,:)=0.
                IF(ierr(i) /= 0) cycle
                DO k=max(kbcon(i),kts+1),ktop(i)
                    dp=100.*(po_cup(i,k+1)-po_cup(i,k))
                    RZenv=0.5*(zuo(i,k+1)+zuo(i,k) - (zdo(i,k+1)+zdo(i,k))*edto(i))
                    S2= cp*tn_cup_x(i,k+1) + g*zo_cup(i,k+1)
                    S1= cp*tn_cup_x(i,k  ) + g*zo_cup(i,k  )
                    Q2= qo_cup_x(i,k+1)
                    Q1= qo_cup_x(i,k  )

                    dqdt(i,k)=         -RZenv*(Q2-Q1)*g/dp
                    dtdt(i,k)= -(1./cp)*RZenv*(S2-S1)*g/dp

                    dqdt(i,k)= dqdt(i,k)+(up_massdetro(i,k)*0.5*(qco (i,k+1)+qco (i,k)-(Q2+Q1)) &
                                + edto(i)*dd_massdetro(i,k)*0.5*(qcdo(i,k+1)+qcdo(i,k)-(Q2+Q1)))*g/dp

                    dtdt(i,k)= dtdt(i,k)+(up_massdetro(i,k)*0.5*(tempco  (i,k+1) + tempco  (i,k)-  &
                                                                (tn_cup_x(i,k+1) + tn_cup_x(i,k))) &
                                + edto(i)*dd_massdetro(i,k)*0.5*(tempcdo (i,k+1) + tempcdo (i,k)-  &
                                                                (tn_cup_x(i,k+1) + tn_cup_x(i,k))) &
                                )*g/dp
                    !print*,"dtdt=",k, dtdt(i,k),zuo(i,k+1),zdo(i,k+1),dqdt(i,k)
                ENDDO
                xk_x(i)=0.
                DO k=max(kbcon(i),kts+1),ktop(i)
                    dp=-(log(100.*po_cup(i,k+1))-log(100.*po_cup(i,k)))      ! no units here
                    xk_x(i)=xk_x(i)+   ( (1.+0.608*qo_x(i,k))*dtdt(i,k) + &
                                                0.608*tn_x(i,k) *dqdt(i,k) )*dp !  units=K m/Pa s2
                    !=> aa3/xk_x will have units of kg/m2/s for the mass flux at cloud base.
                    !print*,"xk_x=",k, xk_x(i),dtdt(i,k),dqdt(i,k)
                ENDDO
            ENDDO
        ENDIF
        !==========================================================
        IF(DICYCLE==5 .or. DICYCLE==2) then

            PHY2:  DO step=fa,fa

                IF(step==bl) then
                    !- version for real cloud-work function by getting the profiles modified only by bl tendencies
                    DO i=its,itf
                        if ( ierr(i) /= 0 ) CYCLE
                        k_free_trop=kpbl(i)
                        !below kbcon -> modify profiles
                        tn_x(i,1:k_free_trop) = tn(i,1:k_free_trop)
                        qo_x(i,1:k_free_trop) = qo(i,1:k_free_trop)
                            !above kbcon -> keep profiles
                        tn_x(i,k_free_trop+1:kte) = t(i,k_free_trop+1:kte)
                        qo_x(i,k_free_trop+1:kte) = q(i,k_free_trop+1:kte)
                    ENDDO
                ELSEIF(step==FA) then
                    !- version for real cloud-work function by getting the profiles modified only by free atmosphere tendencies
                    DO i=its,itf
                        if ( ierr(i) /= 0 ) CYCLE
                        k_free_trop=kpbl(i)
                        !THESE ARE NOT MODIFIED
                        !below kbcon -> keep profiles
                        tn_x(i,1:k_free_trop) = t(i,1:k_free_trop)
                        qo_x(i,1:k_free_trop) = q(i,1:k_free_trop)
                        !THESE MUST BE ONLY FROM ADV + RAD
                        !above kbcon -> modify profiles
                        tn_x(i,k_free_trop+1:kte) = tn(i,k_free_trop+1:kte)-tn_bl(i,k_free_trop+1:kte)+t(i,k_free_trop+1:kte)
                        qo_x(i,k_free_trop+1:kte) = qo(i,k_free_trop+1:kte)-qo_bl(i,k_free_trop+1:kte)+q(i,k_free_trop+1:kte)
                    ENDDO
                ENDIF

                !--- calculate moist static energy, heights, qes, ... only by bl tendencies
                call cup_env(zo,qeso_x,heo_x,heso_x,tn_x,qo_x,po,z1, &
                    psur,ierr,-1,itf,ktf, its,ite, kts,kte)
                !--- environmental values on cloud levels only by bl tendencies
                call cup_env_clev(tn_x,qeso_x,qo_x,heo_x,heso_x,zo,po,qeso_cup_x,qo_cup_x,heo_cup_x,             &
                            us,vs,u_cup,v_cup,                                                   &
                            heso_cup_x,zo_cup,po_cup,gammao_cup_x,tn_cup_x,psur,tsur,  &
                            ierr,z1,itf,ktf,its,ite, kts,kte, host_model)
                do i=its,itf
                    if(ierr(I) /= 0)cycle
                    if(step==bl) then
                        x_add = float(use_excess)*(xlv*zqexec(i)+cp*ztexec(i))
                    else
                        x_add = 0.
                    endif
                    call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup_x(i,kts:kte),hkbo_x(i),k22(i),x_add)
                enddo
                do i=its,itf
                    hco_x (i,:)=0.
                    IF(ierr(i) /= 0)cycle
                    do k=kts,start_level(i)
                        hco_x (i,k) = hkbo_x(i)
                    enddo
                enddo!
                DO i=its,itf
                    if(ierr(i) /= 0)cycle
                    do k=start_level(i)  +1,ktop(i)+1  ! mass cons option
                        denom=(zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
                        if(denom > 0.0)then
                            hco_x (i,k)=(hco_x(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco_x(i,k-1) + &
                                        up_massentro(i,k-1)*heo_x(i,k-1))/ denom
                        else
                            hco_x (i,k)=hco_x(i,k-1)
                        endif
                        hco_x (i,k)=hco_x (i,k) +(1.-p_liq_ice(i,k))*qrco(i,k)*xlf
                    enddo
                    do k=ktop(i)+2,ktf
                        hco_x (i,k)=heso_cup_x(i,k)
                    enddo
                ENDDO
                call get_bouyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
                        ,hco_x,heo_cup_x,heso_cup_x,dbyo_x,zo_cup)

                !--- calculate workfunctions for updrafts
                IF(step==bl) &
                    call cup_up_aa0(aa1_bl,zo_cup,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x &
                                ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

                IF(step==fa) &
                    call cup_up_aa0(aa1_fa,zo_cup,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x &
                                ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)
            ENDDO  PHY2
            DO i=its,itf
                if(ierr(i)/= 0)CYCLE
                aa1_bl(i) = aa1_fa(i)
            ENDDO
        ENDIF
        !---
        !
        !--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
        !
        call cup_dd_edt(ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo, &
            pwo,ccn,pwevo,edtmax,edtmin,maxens2,edtc,psum,psumh, &
            ccnclean,rho,aeroevap,itf,ktf,ipr,jpr,its,ite, kts,kte)

        do 250 iedt=1,maxens2
            do i=its,itf
                if(ierr(i).eq.0)then
                    edto(i)=sigd(i)*edtc(i,iedt)
                    edt (i)=edto(i)
                endif
            enddo
            !
            !--- get the environmental mass flux
            !
            do i=its,itf
                zenv(i,:) = 0.0
                if(ierr(i) /= 0) cycle
                zenv(i,:) = zuo(i,:)-edto(i)*zdo(i,:)
            enddo
            !
            !
            !--- change per unit mass that a model cloud would modify the environment
            !
            !--- 1. in bottom layer
            !
            dellu  =0.
            dellv  =0.
            dellah =0.
            dellat =0.
            dellaq =0.
            dellaqc=0.
            !
            !----------------------------------------------  cloud level ktop
            !
            !- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !
            !----------------------------------------------  cloud level k+2
            !
            !- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
            !
            !----------------------------------------------  cloud level k+1
            !
            !- - - - - - - - - - - - - - - - - - - - - - - - model level k
            !
            !----------------------------------------------  cloud level k
            !
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !      .               .                 .
            !
            !----------------------------------------------  cloud level 3  _cup
            !
            !- - - - - - - - - - - - - - - - - - - - - - - - model level 2
            !
            !----------------------------------------------  cloud level 2  _cup
            !
            !- - - - - - - - - - - - - - - - - - - - - - - - model level 1

            !
            !------------------------------------------------------------------------------------------
            IF(VERT_DISCR == 1) THEN
                do i=its,itf
                    if(ierr(i).eq.0)then
                        do k=kts,ktop(i)
                            ! these three are only used at or near mass detrainment and/or entrainment levels
                            entupk=0.
                            detupk=0.
                            entdoj=0.
                            ! detrainment and entrainment for fowndrafts
                            detdo=edto(i)*dd_massdetro(i,k)
                            entdo=edto(i)*dd_massentro(i,k)
                            ! entrainment/detrainment for updraft
                            entup=up_massentro(i,k)
                            detup=up_massdetro(i,k)
                            ! subsidence by downdrafts only
                            subin=-zdo(i,k+1)*edto(i)
                            subdown=-zdo(i,k)*edto(i)
                            if(k.eq.ktop(i))then
                                detupk=zuo(i,ktop(i))
                                subin=0.
                                subdown=0.
                                detdo=0.
                                entdo=0.
                                entup=0.
                                detup=0.
                            endif
                            totmas=subin-subdown+detup-entup-entdo+ &
                                detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)
                            if(abs(totmas).gt.1.e-6)then
                                write(0,*)'*********************',i,k,totmas
                                !print *,jmin(i),k22(i),kbcon(i),ktop(i)
                                write(0,123)k,subin,subdown,detup,entup, &
                                            detdo,entdo,entupk,detupk
                                123     format(1X,i2,8E12.4)
                                ! call error_fatal ( 'totmas .gt.1.e-6' )
                            endif

                            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                            dellu(i,k) =-(zuo(i,k+1)*(uc (i,k+1)-u_cup(i,k+1) ) -      &
                                        zuo(i,k  )*(uc (i,k  )-u_cup(i,k  ) ) )*g/dp &
                                    +(zdo(i,k+1)*(ucd(i,k+1)-u_cup(i,k+1)) -      &
                                        zdo(i,k  )*(ucd(i,k  )-u_cup(i,k  )) )*g/dp*edto(i)*pgcd

                            dellv(i,k) =-(zuo(i,k+1)*(vc (i,k+1)-v_cup(i,k+1) ) - &
                                        zuo(i,k  )*(vc (i,k  )-v_cup(i,k  ) ) )*g/dp &
                                    +(zdo(i,k+1)*(vcd(i,k+1)-v_cup(i,k+1) ) - &
                                        zdo(i,k  )*(vcd(i,k  )-v_cup(i,k  ) ) )*g/dp*edto(i)*pgcd

                        enddo   ! k

                    endif
                enddo

                do i=its,itf
                    !trash  = 0.0
                    !trash2 = 0.0
                    if(ierr(i).eq.0)then
                        do k=kts,ktop(i)
                            ! these three are only used at or near mass detrainment and/or entrainment levels
                            entupk=0.
                            detupk=0.
                            entdoj=0.
                            ! detrainment and entrainment for downdrafts
                            detdo=edto(i)*dd_massdetro(i,k)
                            entdo=edto(i)*dd_massentro(i,k)
                            ! entrainment/detrainment for updraft
                            entup=up_massentro(i,k)
                            detup=up_massdetro(i,k)
                            ! subsidence by downdrafts only
                            subin=-zdo(i,k+1)*edto(i)
                            subdown=-zdo(i,k)*edto(i)
                            !write(0,*)"down",k,edto(i),detdo,entdo,subin,subdown
                            !if(k.eq.jmin(i))entdoj=edto(i)*zdo(i,k)
                            !if(k.eq.ktop(i))then
                            !   detupk=zuo(i,ktop(i))
                            !   subin=0.; subdown=0.
                            !   detdo=0.; entdo  =0.
                            !   entup=0.; detup  =0.
                            !endif
                            !totmas=subin-subdown+detup-entup-entdo+ &
                            !       detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)

                            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                            dellah(i,k) =-(zuo(i,k+1)*(hco (i,k+1)-heo_cup(i,k+1) ) -                 &
                                        zuo(i,k  )*(hco (i,k  )-heo_cup(i,k  ) ) )*g/dp            &
                                        +(zdo(i,k+1)*(hcdo(i,k+1)-heo_cup(i,k+1) ) -                 &
                                        zdo(i,k  )*(hcdo(i,k  )-heo_cup(i,k  ) ) )*g/dp*edto(i)
                            !---meltglac-------------------------------------------------
                            dellah(i,k) = dellah(i,k) + xlf*((1.-p_liq_ice(i,k))*0.5*(qrco(i,k+1)+qrco(i,k)) &
                                                            - melting(i,k))*g/dp

                            !- check H conservation
                            !trash2 = trash2+ (dellah(i,k))*dp/g

                            !-- take out cloud liquid/ice water for detrainment
                            detup=up_massdetro(i,k)
                            if(.not. USE_C1D  .or. trim(cumulus) == 'mid') then
                                dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                            else
                                if(k == ktop(i)) then
                                    dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                                else
                                    dz=zo_cup(i,k+1)-zo_cup(i,k)
                                    dellaqc(i,k) = zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g
                                endif
                            endif
                            !
                            !---
                            G_rain=  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
                            E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i) ! pwdo < 0 and E_dn must > 0
                            !
                            !write(0,*) "eva=",k,pwdo(i,k),E_dn,zdo(i,k  )
                            !
                            !-- condensation source term = detrained + flux divergence of
                            !-- cloud liquid/ice water (qrco) + converted to rain

                            C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                                zuo(i,k  )* qrco(i,k  )  )*g/dp + G_rain

                            !-- water vapor budget
                            !-- = flux divergence z*(Q_c - Q_env)_up_and_down &
                            !--   - condensation term + evaporation
                            dellaq(i,k) =-(zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) -                 &
                                        zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) )*g/dp            &
                                        +(zdo(i,k+1)*(qcdo(i,k+1)-qo_cup(i,k+1) ) -                 &
                                        zdo(i,k  )*(qcdo(i,k  )-qo_cup(i,k  ) ) )*g/dp*edto(i)    &
                                        - C_up + E_dn

                            !- check water conservation liq+condensed (including rainfall)
                            !trash= trash+ (dellaq(i,k)+dellaqc(i,k)+ G_rain-E_dn)*dp/g

                            !write(3,*)'=>H= ',k,real(trash2,4),real(dellah(i,k),4)
                            !write(4,*)'=>W= ',k,real(trash,4),real(dellaq(i,k),4)
                        enddo   ! k
                        !--- test only with double precision:
                        !write(0,*)'=>H/W-FINAL= ',real(trash2,4),real(trash,4),k22(i),kbcon(i),ktop(i)
                        !if(abs(trash)>1.e-6 .or. abs(trash2) > 1.e-6) then
                        !    write(0,*)'=> not water mass or H cons for deep= ',i,trash,trash2
                        !    !stop 33
                        !endif

                    endif

                enddo
            ELSEIF(VERT_DISCR == 2) THEN
                do i=its,itf
                    if(ierr(i).eq.0)then
                        do k=kts,ktop(i)
                            ! these three are only used at or near mass detrainment and/or entrainment levels
                            entupk=0.
                            detupk=0.
                            entdoj=0.
                            ! detrainment and entrainment for fowndrafts
                            detdo=edto(i)*dd_massdetro(i,k)
                            entdo=edto(i)*dd_massentro(i,k)
                            ! entrainment/detrainment for updraft
                            entup=up_massentro(i,k)
                            detup=up_massdetro(i,k)
                            ! subsidence by downdrafts only
                            subin=-zdo(i,k+1)*edto(i)
                            subdown=-zdo(i,k)*edto(i)
                            if(k.eq.ktop(i))then
                                detupk=zuo(i,ktop(i))
                                subin=0.
                                subdown=0.
                                detdo=0.
                                entdo=0.
                                entup=0.
                                detup=0.
                            endif
                            totmas=subin-subdown+detup-entup-entdo+ &
                                detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)
                            if(abs(totmas).gt.1.e-6)then
                                write(0,*)'*********************',i,k,totmas
                                !print *,jmin(i),k22(i),kbcon(i),ktop(i)
                                write(0,125)k,subin,subdown,detup,entup, &
                                            detdo,entdo,entupk,detupk
                                125     format(1X,i2,8E12.4)
                                ! call error_fatal ( 'totmas .gt.1.e-6' )
                            endif
                            shift = 1
                            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                            dellu(i,k) =-(zuo(i,k+1)*(uc (i,k+1)-us(i,k+shift) ) -      &
                                        zuo(i,k  )*(uc (i,k  )-us(i,max(kts,k-1+shift))) )*g/dp &
                                    +(zdo(i,k+1)*(ucd(i,k+1)-us(i,k+shift)) -      &
                                        zdo(i,k  )*(ucd(i,k  )-us(i,max(kts,k-1+shift))) )*g/dp*edto(i)*pgcd

                            dellv(i,k) =-(zuo(i,k+1)*(vc (i,k+1)-vs(i,k+shift) ) - &
                                        zuo(i,k  )*(vc (i,k  )-vs(i,max(kts,k-1+shift))) )*g/dp &
                                    +(zdo(i,k+1)*(vcd(i,k+1)-vs(i,k+shift) ) - &
                                        zdo(i,k  )*(vcd(i,k  )-vs(i,max(kts,k-1+shift))) )*g/dp*edto(i)*pgcd

                        enddo   ! k

                    endif
                enddo

                do i=its,itf
                    !trash  = 0.0
                    !trash2 = 0.0
                    if(ierr(i).eq.0)then
                        do k=kts,ktop(i)
                            ! these three are only used at or near mass detrainment and/or entrainment levels
                            entupk=0.
                            detupk=0.
                            entdoj=0.
                            ! detrainment and entrainment for downdrafts
                            detdo=edto(i)*dd_massdetro(i,k)
                            entdo=edto(i)*dd_massentro(i,k)
                            ! entrainment/detrainment for updraft
                            entup=up_massentro(i,k)
                            detup=up_massdetro(i,k)
                            ! subsidence by downdrafts only
                            subin=-zdo(i,k+1)*edto(i)
                            subdown=-zdo(i,k)*edto(i)
                            !write(0,*)"down",k,edto(i),detdo,entdo,subin,subdown
                            !if(k.eq.jmin(i))entdoj=edto(i)*zdo(i,k)
                            !if(k.eq.ktop(i))then
                            !   detupk=zuo(i,ktop(i))
                            !   subin=0.; subdown=0.
                            !   detdo=0.; entdo  =0.
                            !   entup=0.; detup  =0.
                            !endif
                            !totmas=subin-subdown+detup-entup-entdo+ &
                            !       detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)

                            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                            !-------------
                            !           dellah(i,k) =-(zuo(i,k+1)*(hco (i,k+1)-heo(i,k+shift) ) -                 &
                            !                          zuo(i,k  )*(hco (i,k  )-heo(i,max(kts,k-1+shift))) )*g/dp            &
                            !                        +(zdo(i,k+1)*(hcdo(i,k+1)-heo(i,k+shift) ) -                 &
                            !                          zdo(i,k  )*(hcdo(i,k  )-heo(i,max(kts,k-1+shift))) )*g/dp*edto(i)
                            !           !---meltglac-------------------------------------------------
                            !           dellah(i,k) = dellah(i,k) + xlf*((1.-p_liq_ice(i,k))*0.5*(qrco(i,k+1)+qrco(i,k)) &
                            !                                              - melting(i,k))*g/dp

                            !--using vert_discr=1 for H
                            dellah(i,k) =-(zuo(i,k+1)*(hco (i,k+1)-heo_cup(i,k+1) ) -                 &
                                        zuo(i,k  )*(hco (i,k  )-heo_cup(i,k  ) ) )*g/dp            &
                                        +(zdo(i,k+1)*(hcdo(i,k+1)-heo_cup(i,k+1) ) -                 &
                                        zdo(i,k  )*(hcdo(i,k  )-heo_cup(i,k  ) ) )*g/dp*edto(i)
                            !---meltglac-------------------------------------------------
                            dellah(i,k) = dellah(i,k) + xlf*((1.-p_liq_ice(i,k))*0.5*(qrco(i,k+1)+qrco(i,k)) &
                                                            - melting(i,k))*g/dp
                            !-------------

                            !- check H conservation
                            !trash2 = trash2+ (dellah(i,k))*dp/g

                            !-- take out cloud liquid/ice water for detrainment
                            detup=up_massdetro(i,k)
                            if(.not. USE_C1D  .or. trim(cumulus) == 'mid') then
                                dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                            else
                                if(k == ktop(i)) then
                                    dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                                else
                                    dz=zo_cup(i,k+1)-zo_cup(i,k)
                                    dellaqc(i,k) = zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g
                                endif
                            endif
                            !
                            !---
                            G_rain=  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
                            E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i) ! pwdo < 0 and E_dn must > 0
                            !
                            !write(0,*) "eva=",k,pwdo(i,k),E_dn,zdo(i,k  )
                            !
                            !-- condensation source term = detrained + flux divergence of
                            !-- cloud liquid/ice water (qrco) + converted to rain

                            C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                                zuo(i,k  )* qrco(i,k  )  )*g/dp + G_rain

                            !-- water vapor budget
                            !-- = flux divergence z*(Q_c - Q_env)_up_and_down &
                            !--   - condensation term + evaporation
                            dellaq(i,k) =-(zuo(i,k+1)*(qco (i,k+1)-qo(i,k+shift) ) -                 &
                                        zuo(i,k  )*(qco (i,k  )-qo(i,max(kts,k-1+shift)) ))*g/dp            &
                                        +(zdo(i,k+1)*(qcdo(i,k+1)-qo(i,k+shift) ) -                 &
                                        zdo(i,k  )*(qcdo(i,k  )-qo(i,max(kts,k-1+shift)) ))*g/dp*edto(i)    &
                                        - C_up + E_dn

                            !- check water conservation liq+condensed (including rainfall)
                            !trash= trash+ (dellaq(i,k)+dellaqc(i,k)+ G_rain-E_dn)*dp/g

                            !write(3,*)'=>H= ',k,real(trash2,4),real(dellah(i,k),4)
                            !write(4,*)'=>W= ',k,real(trash,4),real(dellaq(i,k),4)
                        enddo   ! k
                        !--- test only with double precision:
                        !write(0,*)'=>H/W-FINAL= ',real(trash2,4),real(trash,4),k22(i),kbcon(i),ktop(i)
                        !if(abs(trash)>1.e-6 .or. abs(trash2) > 1.e-6) then
                        !    write(0,*)'=> not water mass or H cons for deep= ',i,trash,trash2
                        !    !stop 33
                        !endif

                    endif

                enddo

            ENDIF ! vertical discretization formulation
            !
            !--- using dellas, calculate changed environmental profiles
            !
            !     if(i.eq.ipr) write(0,*)'DELLAS'
            do k=kts,ktf
                do i=its,itf
                    dellat(i,k)=0.
                    if(ierr(i) /= 0) cycle
                    !
                    XHE(I,K)=(DELLAH(I,K)             )*MBDT(i)+HEO(I,K)
                    XQ (I,K)=(DELLAQ(I,K)+DELLAQC(i,k))*MBDT(i)+QO(I,K)
                    if(XQ(I,K).LE.0.)XQ(I,K)=1.E-08

                    !- do not feed dellat with dellaqc if
                    !- the detrainment of liquid water will be used as
                    !- a source for cloud microphysics
                    if(COUPL_MPHYSICS) then
                        DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xlv*DELLAQ(I,K))
                    else
                        !---meltglac-------------------------------------------------
                        DELLAT (I,K)=(1./cp)*( DELLAH(I,K) - xlv*(DELLAQ(I,K) + DELLAQC(i,k))*(1.+(xlf/xlv)*(1.-p_liq_ice(i,k))))
                        !          DELLAT (I,K)=(1./cp)*( DELLAH(I,K)  -xlv*(DELLAQ(I,K) + DELLAQC(i,k)))

                        !-adding dellaqc to dellaq:
                        DELLAQ (I,K)= DELLAQ(I,K)+DELLAQC(I,K)
                        DELLAQC(I,K)= 0.0
                    endif
                    !---meltglac-------------------------------------------------
                    XT(I,K)=((1./cp)*DELLAH(i,k)-(xlv/cp)*(DELLAQ(i,k)+DELLAQC(i,k)*(1.+(xlf/xlv)*(1.-p_liq_ice(i,k)))))*MBDT(i) &
                    + TN(I,K)
                    !        XT(I,K)=((1./cp)*DELLAH(i,k)-(xlv/cp)*(DELLAQ(i,k)+DELLAQC(i,k)))*MBDT(i)+TN(I,K)

                enddo
            enddo
            do i=its,itf
                if(ierr(i) /= 0) cycle
                !XHKB(I)=(dsubh(i,k22(i))+DELLAH(I,K22(i)))*MBDT+HKBO(I)
                XHE(I,ktf)=HEO(I,ktf)
                XQ(I,ktf)=QO(I,ktf)
                XT(I,ktf)=TN(I,ktf)
                IF(XQ(I,ktf).LE.0.)XQ(I,ktf)=1.E-08
            enddo
            !- new way for defining XHKB
            do i=its,itf
                if(ierr(i) /= 0)CYCLE
                !XHKB(I)= DELLAH(I,K22(i))*MBDT+HKBO(I)
                !-note that HKBO already contains the contribuition from
                !-ztexec and zqexec
                call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),DELLAH (i,kts:kte),DELLAH_aver,k22(i))
                XHKB(I)= DELLAH_aver*MBDT(i) + HKBO(I)
            enddo
            !
            !--- calculate moist static energy, heights, qes
            !
            call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1,psur,ierr,-1,itf,ktf,its,ite,kts,kte)
            !
            !--- environmental values on cloud levels
            !
            call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, xhe_cup,      &
                us,vs,u_cup,v_cup,                                    &
                xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,tsur,   &
                ierr,z1,itf,ktf,its,ite, kts,kte, host_model)
            !
            !--- static control
            !
            !--- moist static energy inside cloud
            !
            do i=its,itf
                xhc(i,:)=0.
                if(ierr(i)/= 0)CYCLE
                do k=kts,start_level(i) !k22(i)
                    xhc(i,k)=xhkb(i)
                enddo
            enddo
            !
            do i=its,itf
                if(ierr(i)/=0)cycle
                do k=start_level(i)  +1,ktop(i)+1  ! mass cons option
                    denom= (xzu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
                    if(denom==0.0) then
                        xhc(i,k)= xhc(i,k-1)
                    else
                        xhc(i,k)=(xhc(i,k-1)*xzu(i,k-1)-.5*up_massdetro(i,k-1)*xhc(i,k-1)+ &
                                                        up_massentro(i,k-1)*xhe(i,k-1)) / denom
                    endif
                    !
                    !- include glaciation effects on XHC
                    !                                   ------ ice content --------
                    xhc (i,k)= xhc (i,k)+ xlf*(1.-p_liq_ice(i,k))*qrco(i,k)
                enddo
                do k=ktop(i)+2,ktf
                    xHC (i,k)=xhes_cup(i,k)
                    xzu (i,k)=0.
                enddo
            enddo
            call get_bouyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
                    ,xhc,xhe_cup,xhes_cup,xdby,xz_cup)
            !
            !--- workfunctions for updraft
            !
            call cup_up_aa0(xaa0,xz_cup,xzu,xdby,GAMMA_CUP,xt_cup, k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

            !if(i.eq.ipr)write(0,*)'xaa0,aa0,aa1 = ',xaa0(1),aa0(1),aa1(1)
            do 200 nens=1,maxens
                do i=its,itf
                if(ierr(i) /= 0) cycle
                !~ xaa0_ens(i,nens)=xaa0(i)

                do k=kts,ktop(i)
                    do nens3=1,maxens3
                        if(nens3.eq.7)then
                        !--- b=0
                            pr_ens(i,nens3)=pr_ens(i,nens3)  &
                                        +pwo(i,k)+edto(i)*pwdo(i,k)
                        !--- b=beta
                        else if(nens3.eq.8)then
                            pr_ens(i,nens3)=pr_ens(i,nens3) &
                                        +pwo(i,k)+edto(i)*pwdo(i,k)
                        !--- b=beta/2
                        else if(nens3.eq.9)then
                            pr_ens(i,nens3)=pr_ens(i,nens3)  &
                                        +pwo(i,k)+edto(i)*pwdo(i,k)

                        else
                            pr_ens(i,nens3)=pr_ens(i,nens3) &
                                        +pwo(i,k)+edto(i)*pwdo(i,k)
                        endif
                    enddo
                enddo
                if(pr_ens(i,7).lt.1.e-6 .and. c0_mid > 0. )then
                    ierr(i)=18
                    ierrc(i)="total normalized condensate too small"
                    do nens3=1,maxens3
                        pr_ens(i,nens3)=0.
                    enddo
                endif
                do nens3=1,maxens3
                    if(pr_ens(i,nens3).lt.1.e-5)then
                        pr_ens(i,nens3)=0.
                    endif
                    enddo
                enddo
            200  continue
            !
            !--- LARGE SCALE FORCING
            !
            do i=its,itf
                ierr2(i)=ierr(i)
                ierr3(i)=ierr(i)
            enddo
            !
            !--- calculate cloud base mass flux
            !
            if(trim(cumulus) == 'deep') &
                call cup_forcing_ens_3d(itf,ktf,its,ite, kts,kte,ens4,ensdim,ichoice,maxens,maxens2,maxens3 &
                                ,ierr,ierr2,ierr3,k22,kbcon,ktop         &
                                ,xland1,aa0,aa1,xaa0,mbdt,dtime          &
                                ,xf_ens,mconv,qo                         &
                                ,po_cup,omeg,zdo,zuo,pr_ens,edto         &
                                ,tau_ecmwf,aa1_bl,xf_dicycle, xk_x)
            !
            if(trim(cumulus) == 'mid') &
                call cup_forcing_ens_3d_mid(aa0,aa1,xaa0,mbdt,dtime,ierr          &
                                    ,po_cup,ktop,k22,kbcon,kpbl,ichoice,maxens,maxens3   &
                                    ,itf,ktf,its,ite, kts,kte                      &
                                    ,tau_ecmwf,aa1_bl,xf_dicycle                   &
                                    ,dhdt,xff_mid,zws,hc,hco,he_cup,heo_cup)

            do k=kts,ktf
                do i=its,itf
                    if(ierr(i) /= 0) cycle
                    pwo_eff(i,k)=pwo(i,k)+edto(i)*pwdo(i,k)
                enddo
            enddo
        250  continue
        !
        !--- Include kinetic energy dissipation converted to heating
        !
        call ke_to_heating(itf,ktf,its,ite, kts,kte,ktop,ierr &
                    ,po_cup,us,vs,dellu,dellv,dellat)
        !
        !--- FEEDBACK
        !
        call cup_output_ens_3d(cumulus,xff_mid,xf_ens,ierr,dellat,dellaq, &
            dellaqc,outt, outq,outqc,zuo,pre,pwo_eff,xmb,ktop,&
            maxens2,maxens,ierr2,ierr3,                  &
            pr_ens,maxens3,ensdim,sig,xland1,                 &
            ichoice,ipr,jpr,itf,ktf,its,ite, kts,kte,         &
            xf_dicycle,outu,outv,dellu,dellv,dtime,po_cup,kbcon )


        !--- get the net precipitation flux (after downdraft evaporation) 
        call get_precip_fluxes(cumulus,klcl,kbcon,ktop,k22,ierr,xland,pre,xmb  &
                        ,pwo,pwavo,edto,pwevo,pwdo,t_cup,tempco          &
                        ,prec_flx,evap_flx                               &
                        ,itf,ktf,its,ite, kts,kte)
        !--- get the total (deep+congestus) evaporation flux for output (units kg/kg/s)
        !
        do i=its,itf
            if(ierr(i) /= 0) cycle
            do k=kts,ktop(i)
                dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                !--- add congestus and deep plumes, and convert to kg/kg/s
                revsu_gf(i,k) = revsu_gf(i,k) + evap_flx(i,k)*g/dp
                prfil_gf(i,k) = prfil_gf(i,k) + prec_flx(i,k)*g/dp
            enddo
        enddo

        !-------------------------- not in use ------------------------------------------------------!
        !--- get cloud fraction
        !
        ! do i=its,itf
        !    clfrac(i,:)=0.
        !    if(ierr(i) /= 0) cycle
        !    dummy1(kts:ktf) = xmb(i)* zuo(i,kts:ktf)
        !    dummy2(kts:ktf) = 100.*po_cup(i,kts:ktf)
        !    call get_cloud_fraction(ktf,kts,ktf						   &
        !     ,dummy2(kts:ktf),zo_cup(i,kts:ktf),tn_cup(i,kts:ktf),qo_cup(i,kts:ktf) &
        !     ,qco (i,kts:ktf),  qrco(i,kts:ktf),  dummy1(kts:ktf),clfrac(i,kts:ktf) )
        ! enddo
        !--------------------------------------------------------------------------------------------!
        !
        !- for tracer convective transport
        do i=its,itf
            if(ierr(i) == 0) then
                do k=kts,ktop(i)+1
                    !clwup5d     (i,k) = qrco         (i,k) !ice/liquid water
                    tup          (i,k) = (1./cp)*(hco(i,k)-g*zo_cup(i,k)-xlv*qco(i,k))!in-updraft temp
                enddo
            endif
        enddo

        !--------------------------------------------------------------------------------------------!
        !- section for atmospheric composition
        !--------------------------------------------------------------------------------------------!
        IF(USE_TRACER_TRANSP==1)  THEN

            !-0) convert mass fluxes, etc...
            do i=its,itf
                if(ierr(i) /= 0) cycle
                pwavo       (i)   = xmb(i)*pwavo       (i)
                pwevo       (i)   = xmb(i)*pwevo       (i)
                zuo         (i,:) = xmb(i)*zuo         (i,:)
                zdo         (i,:) = xmb(i)*zdo         (i,:)
                pwo         (i,:) = xmb(i)*pwo         (i,:)
                pwdo        (i,:) = xmb(i)*pwdo        (i,:)
                up_massentro(i,:) = xmb(i)*up_massentro(i,:)
                up_massdetro(i,:) = xmb(i)*up_massdetro(i,:)
                dd_massentro(i,:) = xmb(i)*dd_massentro(i,:)
                dd_massdetro(i,:) = xmb(i)*dd_massdetro(i,:) 
                zenv        (i,:) = xmb(i)*zenv        (i,:)
            enddo

            !-1) get mass mixing ratios at the cloud levels

            call cup_env_clev_chem(mtp,se_chem,se_cup_chem,ierr,itf,ktf,its,ite, kts,kte)

            !-2) determine in-cloud tracer mixing ratios
            !
            ! a) chem - updraft
            !- note: here "sc_up_chem" stores the total in-cloud tracer mixing ratio (i.e., including the portion
            !        embedded in the condensates).
            call get_incloud_sc_chem_up(cumulus,mtp,se_chem,se_cup_chem,sc_up_chem,pw_up_chem,tot_pw_up_chem      &
                            ,zo_cup,rho,po,po_cup,qrco,tempco,pwo,zuo,up_massentro,up_massdetro                &
                            ,vvel2d,vvel1d,start_level,k22,kbcon,ktop,klcl,ierr,xland,itf,ktf,its,ite, kts,kte)

            ! b) chem - downdraft
            call get_incloud_sc_chem_dd(cumulus,mtp,se_chem,se_cup_chem,sc_dn_chem,pw_dn_chem ,pw_up_chem,sc_up_chem &
                            ,tot_pw_up_chem,tot_pw_dn_chem                                                       & 
                            ,zo_cup,rho,po_cup,qrcdo,pwdo,pwevo,edto,zdo,dd_massentro,dd_massdetro ,pwavo,pwo     &
                            ,jmin,ierr,itf,ktf,its,ite, kts,kte)
            !
            !-3) determine the vertical transport including mixing, scavenging and evaporation
            !
            !---a) change per unit mass that a model cloud would modify the environment
            do i=its,itf
                if(ierr(i) /= 0) cycle

                !- flux form + source/sink terms + time explicit + FCT
                IF(USE_FLUX_FORM == 1 .and. alp1 == 0. ) THEN 

                    if(use_fct == 0 ) then 
                        do k=kts,ktop(i)
                            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                            out_chem(:,i,k) =-(zuo(i,k+1)*(sc_up_chem(:,i,k+1)-se_cup_chem(:,i,k+1) ) -                 &
                                                zuo(i,k  )*(sc_up_chem(:,i,k  )-se_cup_chem(:,i,k  ) ))*g/dp             &
                                                +(zdo(i,k+1)*(sc_dn_chem(:,i,k+1)-se_cup_chem(:,i,k+1) ) -                 &
                                                zdo(i,k  )*(sc_dn_chem(:,i,k  )-se_cup_chem(:,i,k  ) ))*g/dp*edto(i)
                        enddo

                    else 

                        !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
                        sub_tend  = 0.
                        trcflx_in = 0. 
                        dtime_max = dtime
                        massflx  (i,kts)=0.

                        do k=kts+1,ktop(i)+1
                            dp  	       = 100.*(po_cup(i,k)-po_cup(i,k+1))
                            trcflx_in (:,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))*se_cup_chem(:,i,k) !* xmb(i)
                            massflx   (i,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))		       !* xmb(i)
                            dtime_max=min(dtime_max,.5*dp) 
                        enddo
                        !- if dtime_max<dtime => needs a loop to update from t to t+dtime (check this!)
                        !if( dtime_max < dtime ) stop "dtime_max < dtime in GF scheme"

                        do ispc=1,mtp
                            call fct1d3 (ktop(i),kte,dtime_max,po_cup(i,:),se_chem(ispc,i,:),massflx(i,:),trcflx_in(ispc,:),sub_tend(ispc,:))
                        enddo

                        do k=kts,ktop(i)
                            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                            out_chem(:,i,k) = -(zuo(i,k+1)*(sc_up_chem(:,i,k+1)) - zuo(i,k)*(sc_up_chem(:,i,k)))*g/dp  	&
                                            +(zdo(i,k+1)*(sc_dn_chem(:,i,k+1)) - zdo(i,k)*(sc_dn_chem(:,i,k)))*g/dp*edto(i)

                            !- update with the subsidence term from FCT scheme
                            out_chem(:,i,k) = out_chem(:,i,k) + sub_tend(:,k)

                        enddo
                    endif

                    !- include evaporation
                    if(USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
                        do k=kts,ktop(i)
                            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                            out_chem(:,i,k) = out_chem(:,i,k)    &
                                                - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1)*pw_dn_chem(:,i,k+1))*g/dp !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
                        enddo
                    endif

                    !- include scavenging
                    if(USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
                        do k=kts,ktop(i)
                            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                            out_chem(:,i,k) = out_chem(:,i,k) &
                                            - 0.5*(zuo(i,k)*pw_up_chem(:,i,k)+zuo(i,k+1)*pw_up_chem(:,i,k+1))*g/dp  ! incorporated in rainfall (<0)
                        enddo
                    endif
                ENDIF

                !- flux form + source/sink terms + time explicit/implicit + upstream
                IF(USE_FLUX_FORM == 1 .and. alp1 > 0. ) THEN 

                    alp0=1.-alp1
                    do k=kts,ktop(i)+1
                        fp(k) = 0.5*(zenv(i,k)+abs(zenv(i,k)))
                        fm(k) = 0.5*(zenv(i,k)-abs(zenv(i,k)))
                    enddo

                    do k=kts,ktop(i)
                        dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                        beta1 = dtime*g/dp
                        aa(k) =    alp1*beta1*fm(k)
                        bb(k) = 1.+alp1*beta1*(fp(k)-fm(k+1))
                        cc(k) =   -alp1*beta1*fp(k+1)

                        ddtr(:,k) = se_chem(:,i,k) - (zuo(i,k+1)*sc_up_chem(:,i,k+1) - zuo(i,k)*sc_up_chem(:,i,k))*beta1      &
                                                    + (zdo(i,k+1)*sc_dn_chem(:,i,k+1) - zdo(i,k)*sc_dn_chem(:,i,k))*beta1*edto(i)

                        !- include evaporation
                        if(USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
                            out_chem(:,i,k) = out_chem(:,i,k)    &
                                            - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1)*pw_dn_chem(:,i,k+1))*beta1 !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
                        endif

                        !- include scavenging
                        if(USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
                            out_chem(:,i,k) = out_chem(:,i,k) &
                                            - 0.5*(zuo(i,k)*pw_up_chem(:,i,k)+zuo(i,k+1)*pw_up_chem(:,i,k+1))*beta1  ! incorporated in rainfall (<0)
                        endif

                        ddtr(:,k) = ddtr(:,k) + out_chem(:,i,k) + & 
                                alp0*beta1*(-fm(k)*se_chem(:,i,max(kts,k-1)) +(fm(k+1)-fp(k))*se_chem(:,i,k) +fp(k+1)*se_chem(:,i,k+1))

                    enddo
                    do ispc = 1, mtp
                        call tridiag (ktop(i),aa (kts:ktop(i)),bb (kts:ktop(i)),cc (kts:ktop(i)),ddtr (ispc,kts:ktop(i)))
                        out_chem(ispc,i,kts:ktop(i))=(ddtr(ispc,kts:ktop(i))-se_chem(ispc,i,kts:ktop(i)))/dtime
                    enddo	    
                ENDIF

                !- flux form + source/sink terms + time explicit + upstream with anti-diffusion step (Smolarkiewicz 1983)
                IF(USE_FLUX_FORM == 2 .or. USE_FLUX_FORM == 3) THEN 
                    if(USE_FLUX_FORM == 2)  lstep = -1 ! upstream + anti-diffusion step
                    if(USE_FLUX_FORM == 3)  lstep =  1 ! only upstream
                    alp0=1.

                    if(ierr(i) /= 0) cycle
                    !--- Zenv here have the following reference:  < 0 => downward motion
                    zenv(i,:) = -(zuo(i,:)-edto(i)*zdo(i,:))

                    do istep=1,lstep,-2

                        if(istep == 1) then 
                            ddtr_upd(:,:) = se_chem(:,i,:)
                            do k=kts,ktop(i)+1
                                fp_mtp(:,k) = 0.5*(zenv(i,k)+abs(zenv(i,k)))
                                fm_mtp(:,k) = 0.5*(zenv(i,k)-abs(zenv(i,k)))
                            enddo
                        else
                            ddtr_upd(:,kts:ktop(i)+1) = se_chem(:,i,kts:ktop(i)+1) + out_chem(:,i,kts:ktop(i)+1)*dtime
                            zenv_diff(:,kts) = 0.
                            do k=kts,ktop(i)+1
                                dz =  zo_cup(i,k+1)-zo_cup(i,k) 
                                zenv_diff (:,k+1) = 1.08*( dz*abs(zenv(i,k+1)) - dtime * zenv(i,k+1)**2 ) &
                                            * (ddtr_upd(:,k+1) - ddtr_upd(:,k))               &
                                                    /((ddtr_upd(:,k+1) + ddtr_upd(:,k) + 1.e-16)*dz) 
                            enddo
                            do k=kts,ktop(i)+1
                                fp_mtp(:,k) = 0.5*(zenv_diff(:,k)+abs(zenv_diff(:,k)))
                                fm_mtp(:,k) = 0.5*(zenv_diff(:,k)-abs(zenv_diff(:,k)))
                            enddo
                        endif
               
                        do k=kts,ktop(i)
                            dp=-100.*(po_cup(i,k)-po_cup(i,k+1))
                            beta1 = dtime*g/dp
                            ddtr(:,k) = ddtr_upd(:,k) + alp0*beta1*( &
                                                    (fp_mtp(:,k+1)*ddtr_upd(:,k)           +fm_mtp(:,k+1)*ddtr_upd(:,k+1)) &
                                                    -(fp_mtp(:,k  )*ddtr_upd(:,max(kts,k-1))+fm_mtp(:,k  )*ddtr_upd(:,k  )) )
                        enddo
                        do ispc = 1, mtp
                            out_chem(ispc,i,kts:ktop(i))=(ddtr(ispc,kts:ktop(i))-se_chem(ispc,i,kts:ktop(i)))/dtime
                        enddo

                    enddo ! anti-diff steps

                    do k=kts,ktop(i)
                        dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                        beta1 = g/dp

                        out_chem(:,i,k)  =   out_chem(:,i,k)                                                           &
                                        - (zuo(i,k+1)*sc_up_chem(:,i,k+1) - zuo(i,k)*sc_up_chem(:,i,k))*beta1      &
                                            + (zdo(i,k+1)*sc_dn_chem(:,i,k+1) - zdo(i,k)*sc_dn_chem(:,i,k))*beta1*edto(i)

                        !- include evaporation
                        if(USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
                            out_chem(:,i,k) = out_chem(:,i,k)    &
                                            - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1)*pw_dn_chem(:,i,k+1))*beta1 !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
                        endif

                        !- include scavenging
                        if(USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
                            out_chem(:,i,k) = out_chem(:,i,k) &
                                            - 0.5*(zuo(i,k)*pw_up_chem(:,i,k)+zuo(i,k+1)*pw_up_chem(:,i,k+1))*beta1  ! incorporated in rainfall (<0)
                        endif

                    enddo

                ENDIF

                !--- check mass conservation for tracers 
                do ispc = 1, mtp
                    evap_  (:) = 0.
                    wetdep_(:) = 0.
                    residu_(:) = 0.
                    do k=kts,ktop(i)
                        dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                        evap   =  -0.5*(zdo(i,k)*pw_dn_chem(ispc,i,k)+zdo(i,k+1)*pw_dn_chem(ispc,i,k+1))*g/dp*edto(i) 
                        wetdep =   0.5*(zuo(i,k)*pw_up_chem(ispc,i,k)+zuo(i,k+1)*pw_up_chem(ispc,i,k+1))*g/dp

                        evap_  (ispc) =   evap_  (ispc) + evap  *dp/g 
                        wetdep_(ispc) =   wetdep_(ispc) + wetdep*dp/g 
                        residu_(ispc) =   residu_(ispc) + (wetdep - evap)*dp/g
                        
                    enddo
                    if(residu_(ispc) < 0.) then 
                        beta1 = g/(po_cup(i,kts)-po_cup(i,ktop(i)+1)) 
                        do k=kts,ktop(i)
                            out_chem(ispc,i,k)=out_chem(ispc,i,k)+residu_(ispc)*beta1
                        enddo
                    endif
                enddo


            enddo ! loop 'i'

            !-4) convert back mass fluxes, etc...
            do i=its,itf
                if(ierr(i) /= 0) cycle
                pwavo       (i)   = pwavo (i)  / (xmb(i) + 1.e-16)
                pwevo       (i)   = pwevo (i)  / (xmb(i) + 1.e-16)
                zuo         (i,:) = zuo   (i,:)/ (xmb(i) + 1.e-16)
                zdo         (i,:) = zdo   (i,:)/ (xmb(i) + 1.e-16)
                pwo         (i,:) = pwo   (i,:)/ (xmb(i) + 1.e-16)
                pwdo        (i,:) = pwdo  (i,:)/ (xmb(i) + 1.e-16)
                up_massentro(i,:) = up_massentro(i,:)/ (xmb(i) + 1.e-16)
                up_massdetro(i,:) = up_massdetro(i,:)/ (xmb(i) + 1.e-16)
                dd_massentro(i,:) = dd_massentro(i,:)/ (xmb(i) + 1.e-16)
                dd_massdetro(i,:) = dd_massdetro(i,:)/ (xmb(i) + 1.e-16)
                zenv        (i,:) = zenv  (i,:)/ (xmb(i) + 1.e-16)
            enddo

            !--------------------------------------------------------------------------------------------!
        ENDIF !- end of section for atmospheric composition
        !--------------------------------------------------------------------------------------------!

        !- for debug/diag
        if(trim(cumulus) == 'deep') then
            do i=its,itf
                !if(ierr(i) /= 0) cycle
                aa0_    (i)  = aa0(i)
                aa1_    (i)  = edto(i) !aa1(i)
                aa2_    (i)  = aa2(i)
                aa3_    (i)  = aa3(i)
                aa1_bl_ (i)  = aa1_bl(i)
                aa1_cin_(i)  = cin1(i)
                tau_bl_ (i)  = tau_bl   (i)
                tau_ec_ (i)  = tau_ecmwf(i)
            ENDDO
        endif
        !
        !- begin: for GATE soundings-------------------------------------------
        if(wrtgrads) then
            if(trim(cumulus) == 'deep') then ; cty='1'; nvarbegin =  1 ;endif
            if(trim(cumulus) == 'mid' ) then ; cty='3'; nvarbegin = 51 ;endif
            do i=its,itf
                !if(ierr(i).eq.0) then
                !- 2-d section
                do k=kts,max(1,ktop(i)+1)

                    dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                    E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i)*86400.*xlv/cp*xmb(i) ! pwdo < 0 and E_dn must > 0
                    C_up  = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) - zuo(i,k  )* qrco(i,k  )  )*g/dp &
                            +0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
                    C_up = - C_up*86400.*xlv/cp*xmb(i)

                    trash =-(zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) -                 &
                    zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) )*g/dp            
                    trash2=+(zdo(i,k+1)*(qcdo(i,k+1)-qo_cup(i,k+1) ) -                 &
                    zdo(i,k  )*(qcdo(i,k  )-qo_cup(i,k  ) ) )*g/dp*edto(i)    
                    
                    trash  = trash *86400.*xlv/cp*xmb(i); trash2 = trash2*86400.*xlv/cp*xmb(i)

                    nvar=nvarbegin
                    ! call set_grads_var(jl,k,nvar,out_chem(1,i,k)*86400,"outchem"//cty ,' outchem','3d')
                    ! call set_grads_var(jl,k,nvar,sc_up_chem(1,i,k),"scup"//cty ,' sc_chem','3d')
                    ! call set_grads_var(jl,k,nvar,sc_dn_chem(1,i,k),"scdn"//cty ,' sc_chem','3d')
                    ! call set_grads_var(jl,k,nvar,massi,"mi"//cty ,' initial mass','2d')
                    ! call set_grads_var(jl,k,nvar,massf,"mf"//cty ,' final mass','2d')
                    ! call set_grads_var(jl,k,nvar,se_chem(1,i,k),"se"//cty ,' se_chem','3d')
                    ! call set_grads_var(jl,k,nvar,se_cup_chem(1,i,k),"secup"//cty ,' se_cup_chem','3d')

                    ! call set_grads_var(jl,k,nvar,zuo(i,k),"zup"//cty,'norm m flux up ','3d')
                    ! call set_grads_var(jl,k,nvar,zdo(i,k),"zdn"//cty,'norm m flux dn ','3d')
                    ! call set_grads_var(jl,k,nvar,xmb(i),"xmb"//cty,'m flux up (kg/s/m^2)','2d')
                    ! call set_grads_var(jl,k,nvar,-edto(i)*xmb(i)*zdo(i,k),"mdn"//cty ,'m flux down (kg/s/m^2)','3d')
                    ! call set_grads_var(jl,k,nvar,dellah(i,k)*86400./cp,"delh"//cty ,'dellah K/day','3d')
                    ! call set_grads_var(jl,k,nvar,dellaq(i,k)*86400.*xlv/cp, "dellq"//cty ,'dellaq K/day','3d')
                    ! call set_grads_var(jl,k,nvar,dellaqc(i,k)*86400.*xlv/cp,"dellqc"//cty ,'dellaqc K/day','3d')
                    ! call set_grads_var(jl,k,nvar,xmb(i)*up_massentro(i,k),"upent"//cty ,'up_massentr(kg/s/m^2)','3d')
                    ! call set_grads_var(jl,k,nvar,xmb(i)*up_massdetro(i,k),"updet"//cty ,'up_massdetr(kg/s/m^2)','3d')
                    ! call set_grads_var(jl,k,nvar,outt(i,k)*86400.,"outt"//cty ,'outt K/day','3d')
                    ! call set_grads_var(jl,k,nvar,outq(i,k)*86400.*xlv/cp,"outq"//cty ,'outq K/s','3d')
                    ! call set_grads_var(jl,k,nvar,outqc(i,k)*86400.*xlv/cp,"outqc"//cty ,'outqc K/day','3d')
                    ! call set_grads_var(jl,k,nvar,pre(i)*3600.,"precip"//cty ,'precip mm','2d')
                    ! call set_grads_var(jl,k,nvar,outu(i,k)*86400.,"outu"//cty ,'out_U m/s/day','3d')
                    ! call set_grads_var(jl,k,nvar,outv(i,k)*86400.,"outv"//cty ,'out_V m/s/day','3d')
                    ! call set_grads_var(jl,k,nvar,xmb(i),"xmb"//cty ,'xmb kg/m2/s','2d')
                    ! call set_grads_var(jl,k,nvar,float(ierr(i)),"ierr"//cty ,'ierr #','2d')
                    ! call set_grads_var(jl,k,nvar,xmb(i)*dd_massentro(i,k),"ddent"//cty ,'dd_massentr(kg/s/m^2)','3d')
                    ! call set_grads_var(jl,k,nvar,xmb(i)*dd_massdetro(i,k),"dddet"//cty ,'dd_massdetr(kg/s/m^2)','3d')
                    ! !call set_grads_var(jl,k,nvar,QCUP(i,k),"qcup"//cty ,'C_UP','3d')
                    ! call set_grads_var(jl,k,nvar,t_cup(i,k)-273.15,"te"//cty ,' K','3d')
                    ! call set_grads_var(jl,k,nvar,1000.*q_cup(i,k),"qe"//cty ,' kg kg-1','3d')
                    ! call set_grads_var(jl,k,nvar,he_cup(i,k),"he"//cty ,' he','3d')
                    ! call set_grads_var(jl,k,nvar,HKB(i),"hkb"//cty ,' H','2d')
                    ! call set_grads_var(jl,k,nvar,HKB(i),"hkb"//cty ,' H','2d')
                    ! call set_grads_var(jl,k,nvar,z_cup(i,k22(i)),"zs"//cty ,' m','2d')
                    ! call set_grads_var(jl,k,nvar,z_cup(i,kbcon(i)),"zbcon"//cty ,' m','2d')
                    ! call set_grads_var(jl,k,nvar,z_cup(i,ktop(i)),"ztop"//cty ,' m','2d')
                    ! call set_grads_var(jl,k,nvar,z_cup(i,klcl(i)),"zlcl"//cty ,' m','2d')
                    ! call set_grads_var(jl,k,nvar,zws(i),"ws"//cty ,' m/s','2d')
                    ! call set_grads_var(jl,k,nvar,clfrac(i,k),"clfrac"//cty ,'shcf #','3d')
                    ! call set_grads_var(jl,k,nvar,vvel2d(i,k),"W2d"//cty ,'W /m/s','3d')
                    ! call set_grads_var(jl,k,nvar,vvel1d(i),"W1d"//cty ,'W1s /m/s','2d')
                    ! call set_grads_var(jl,k,nvar,entr_rate_2d(i,k),"entr"//cty ,' m-1','3d')
                    ! call set_grads_var(jl,k,nvar,cd(i,k),"detr"//cty ,' m-1','3d')
                    ! call set_grads_var(jl,k,nvar,pwdo(i,k),"pwd"//cty ,' xx','3d')
                    ! call set_grads_var(jl,k,nvar,edto(i),"edt"//cty ,'edt kg/m2/s','2d')
                    ! call set_grads_var(jl,k,nvar,E_DN,"EVAP"//cty ,' xx','3d')
                    ! call set_grads_var(jl,k,nvar,C_UP,"CUP"//cty ,' xx','3d')
                    ! call set_grads_var(jl,k,nvar,trash,"TUP"//cty ,' xx','3d')
                    ! call set_grads_var(jl,k,nvar,trash2,"TDN"//cty ,' xx','3d')

                    ! call set_grads_var(jl,k,nvar,p_liq_ice(i,k),"pli"//cty ,'#','3d')
                    ! call set_grads_var(jl,k,nvar,melting_layer(i,k),"cpli"//cty ,'#','3d')
                    ! call set_grads_var(jl,k,nvar,t(i,k),"t"//cty ,'temp K','3d')
                    ! call set_grads_var(jl,k,nvar,tn(i,k),"tn"//cty ,'temp K','3d')
                    ! call set_grads_var(jl,k,nvar,1000.*q(i,k),"q"//cty ,'q g/kg','3d')
                    ! call set_grads_var(jl,k,nvar,1000.*qo(i,k),"qn"//cty ,'q g/kg','3d')
                    ! call set_grads_var(jl,k,nvar,1000.*qrco(i,k),"qrc"//cty ,'q g/kg','3d')
                    ! call set_grads_var(jl,k,nvar,1000.*(q(i,k)+outq(i,k)*dtime),"qnc"//cty ,'q upd conv g/kg','3d')
                    ! call set_grads_var(jl,k,nvar,1000.*(qo(i,k)+outq(i,k)*dtime),"qnall"//cty ,'q upd all g/kg','3d')

                    !~ call set_grads_var(jl,k,nvar,aa0(i),"a0"//cty,'aa0','2d')
                    !~ call set_grads_var(jl,k,nvar,aa1_fa(i),"aa1fa"//cty,'aa1fa','2d')
                    !~ call set_grads_var(jl,k,nvar,aa1_bl(i),"aa1bl"//cty,'aa1bl','2d')
                    !~ call set_grads_var(jl,k,nvar,aa0_bl(i),"aa0bl"//cty,'aa0bl','2d')
                    !~ call set_grads_var(jl,k,nvar,aa1(i),"a1"//cty,'aa1','2d')
                    !~ call set_grads_var(jl,k,nvar,aa1(i)/(1.e-6+tau_ecmwf(i)),"mb13"//cty,'aa0','2d')
                    !~ call set_grads_var(jl,k,nvar,xaa0(i),"xa0"//cty,'xaa0','2d')
                    !~ call set_grads_var(jl,k,nvar,(XAA0(I)-AA1(I))/MBDT(I),"xk"//cty,'xk','2d')

                enddo
                ! if(wrtgrads) then
                !     call wrt_bin_ctl(1,kte,po(1,1:kte),cumulus)
                ! endif
            enddo
        endif
        !- end  : for GATE soundings-------------------------------------------

    END SUBROUTINE CUP_gf_GEOS5

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