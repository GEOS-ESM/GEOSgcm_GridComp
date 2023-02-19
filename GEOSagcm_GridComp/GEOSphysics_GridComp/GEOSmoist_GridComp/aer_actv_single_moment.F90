MODULE Aer_Actv_Single_Moment
!
#include "MAPL_Generic.h"

      USE ESMF
      USE MAPL
      USE aer_cloud, only: AerProps
!-------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      PUBLIC ::  Aer_Activation, USE_BERGERON, USE_AEROSOL_NN, R_AIR
      PRIVATE

       ! Real kind for activation.
       integer,public,parameter :: AER_R4 = SELECTED_REAL_KIND(6,37)
       integer,public,parameter :: AER_R8 = SELECTED_REAL_KIND(15,307)
       integer,public,parameter :: AER_PR = AER_R8

       real        , parameter :: R_AIR     =  3.47e-3 !m3 Pa kg-1K-1
       real(AER_PR), parameter :: zero_par  =  1.e-6   ! small non-zero value
       real(AER_PR), parameter :: ai        =  0.0000594D0
       real(AER_PR), parameter :: bi        =  3.33D0
       real(AER_PR), parameter :: ci        =  0.0264D0
       real(AER_PR), parameter :: di        =  0.0033D0

       real(AER_PR), parameter :: betai     = -2.262D+3
       real(AER_PR), parameter :: gamai     =  5.113D+6
       real(AER_PR), parameter :: deltai    =  2.809D+3
       real(AER_PR), parameter :: densic    =  917.D0   !Ice crystal density in kgm-3

       real, parameter :: NN_LAND     = 300.0e6
       real, parameter :: NN_OCEAN    = 100.0e6
       real, parameter :: NN_MIN      =  30.0e6

       LOGICAL  :: USE_BERGERON, USE_AEROSOL_NN
      CONTAINS          

!>----------------------------------------------------------------------------------------------------------------------
!>----------------------------------------------------------------------------------------------------------------------

      SUBROUTINE Aer_Activation(IM,JM,LM, q, t, plo, ple, zlo, zle, qlcn, qicn, qlls, qils, &
                                       sh, evap, kpbl, omega, FRLAND, USE_AERO_BUFFER, &
                                       AeroProps, aero_aci, NACTL, NACTI, NWFA)
      IMPLICIT NONE
      integer, intent(in)::IM,JM,LM
      TYPE(AerProps), dimension (IM,JM,LM),intent(inout)  :: AeroProps
      type(ESMF_State)            ,intent(inout) :: aero_aci
      real, dimension (IM,JM,LM)  ,intent(in ) :: plo ! Pa
      real, dimension (IM,JM,0:LM),intent(in ) :: ple ! Pa
      real, dimension (IM,JM,LM)  ,intent(in ) :: q,t,omega,zlo, qlcn, qicn, qlls, qils
      real, dimension (IM,JM,0:LM),intent(in ) :: zle
      real, dimension (IM,JM)     ,intent(in ) :: FRLAND
      real, dimension (IM,JM)     ,intent(in ) :: sh, evap, kpbl     
      logical                     ,intent(in ) :: USE_AERO_BUFFER
      
 
      real, dimension (IM,JM,LM),intent(OUT) :: NACTL,NACTI, NWFA
      
      real(AER_PR), allocatable, dimension (:) :: sig0,rg,ni,bibar,nact 
      real(AER_PR), dimension (IM,JM)  :: naer_cb, zws
      real(AER_PR)                     :: wupdraft,tk,press,air_den,QC,QL,WC,BB,RAUX

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
      REAL                            :: aux1,aux2,aux3,hfs,hfl, nfaux

      integer :: n_modes
      REAL :: numbinit
      integer :: i,j,k,n,rc

      character(len=ESMF_MAXSTR)              :: IAm="Aer_Activation"
      integer                                 :: STATUS

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

          call ESMF_AttributeGet(aero_aci, name='number_of_aerosol_modes', value=n_modes, __RC__)

          if (n_modes > 0) then

              allocate( sig0(n_modes), __STAT__)
              allocate(   rg(n_modes), __STAT__)
              allocate(   ni(n_modes), __STAT__)
              allocate(bibar(n_modes), __STAT__)
              allocate( nact(n_modes), __STAT__)

              allocate(aero_aci_modes(n_modes), __STAT__)
              call ESMF_AttributeGet(aero_aci, name='aerosol_modes', itemcount=n_modes, valuelist=aero_aci_modes, __RC__)

              call ESMF_AttributeGet(aero_aci, name='air_pressure_for_aerosol_optics', value=aci_field_name, __RC__)
              if (aci_field_name /= '') then
                  call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
                  aci_ptr_3d = PLE
              end if

              call ESMF_AttributeGet(aero_aci, name='air_temperature', value=aci_field_name, __RC__)
              if (aci_field_name /= '') then
                  call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
                  aci_ptr_3d = T
              end if

              call ESMF_AttributeGet(aero_aci, name='fraction_of_land_type', value=aci_field_name, __RC__)
              if (aci_field_name /= '') then
                  call MAPL_GetPointer(aero_aci, aci_ptr_2d, trim(aci_field_name), __RC__)
                  aci_ptr_2d = FRLAND
              end if
     
              if (USE_AERO_BUFFER) then
                 allocate(buffer(im,jm,lm,n_modes,8), __STAT__)
              end if

              ACTIVATION_PROPERTIES: do n = 1, n_modes
                 call ESMF_AttributeSet(aero_aci, name='aerosol_mode', value=trim(aero_aci_modes(n)), __RC__)
                 
                 ! execute the aerosol activation properties method 
                 call ESMF_MethodExecute(aero_aci, label='aerosol_activation_properties', userRC=ACI_STATUS, RC=STATUS)
                 VERIFY_(ACI_STATUS)
                 VERIFY_(STATUS)

                 ! copy out aerosol activation properties
                 call ESMF_AttributeGet(aero_aci, name='aerosol_number_concentration', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_num, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='aerosol_dry_size', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_dgn, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='width_of_aerosol_mode', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_sigma, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='aerosol_density', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_density, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='aerosol_hygroscopicity', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_hygroscopicity, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='fraction_of_dust_aerosol', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_f_dust, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='fraction_of_soot_aerosol', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_f_soot, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='fraction_of_organic_aerosol', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_f_organic, trim(aci_field_name), __RC__)

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

              if (USE_AERO_BUFFER) then
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

                 deallocate(buffer, __STAT__)
              end if
              
             
              do k = 1, LM
               do j = 1, JM
                do i = 1, IM
                nfaux =  0.0
                 do n = 1, n_modes
                   if (AeroProps(i,j,k)%kap(n) .gt. 0.4) then 
                            nfaux =  nfaux + AeroProps(i,j,k)%num(n)
                   end if                           
                 end do !modes
                 NWFA(I, J, K)  =  nfaux
                end do
               end do 
              end do 
             

              deallocate(aero_aci_modes, __STAT__)

!----- aerosol activation (single-moment uphysics)      
      do j = 1, JM
         do i = 1, IM
            aux1=  PLE(i,j,LM)/(287.04*(T(i,j,LM)*(1.+0.608*Q(i,j,LM)))) ! air_dens (kg m^-3)
            hfs = -SH  (i,j) ! W m^-2
            hfl = -EVAP(i,j) ! kg m^-2 s^-1
            aux2= (hfs/MAPL_CP + 0.608*T(i,j,LM)*hfl)/aux1 ! buoyancy flux (h+le)
            aux3=  ZLE(i,j,kpbli(i,j))           ! pbl height (m)
            !-convective velocity scale W* (m/s)
            ZWS(i,j) = max(0.,0.001-1.5*0.41*MAPL_GRAV*aux2*aux3/T(i,j,LM))
            ZWS(i,j) = 1.2*ZWS(i,j)**0.3333 ! m/s           
      enddo; enddo

      !--- activated aerosol # concentration for liq/ice phases (units: m^-3)
      numbinit = 0.
      WC       = 0.
      BB       = 0.
      RAUX     = 0.
           
      !--- determing aerosol number concentration at cloud base
      DO j=1,JM
        Do i=1,IM 
             k            = kpbli(i,j)
             naer_cb(i,j) = zero_par
             tk           = T(i,j,k)              ! K
             press        = plo(i,j,k)            ! Pa     
             air_den      = press*28.8e-3/8.31/tk ! kg/m3
             DO n=1,n_modes
                if (AeroProps(i,j,k)%dpg(n) .ge. 0.5e-6) &
                naer_cb(i,j)= naer_cb(i,j) + AeroProps(i,j,k)%num(n)         
             ENDDO
             naer_cb(i,j)= naer_cb(i,j) * air_den * 1.e+6  ! #/cm3
             naer_cb(i,j)= max(0.1,min(naer_cb(i,j),100.0))
      ENDDO;ENDDO
     
      DO k=LM,1,-1
       DO j=1,JM
        Do i=1,IM
              
              tk                 = T(i,j,k)                         ! K
              press              = plo(i,j,k)                       ! Pa   
              air_den            = press*28.8e-3/8.31/tk            ! kg/m3
              qc                 = (qicn(i,j,k)+qils(i,j,k))*1.e+3  ! g/kg
              ql                 = (qlcn(i,j,k)+qlls(i,j,k))*1.e+3  ! g/kg
              
              IF( plo(i,j,k) > 34000.0) THEN 
                
                wupdraft           = -9.81*air_den*omega(i,j,k)     ! m/s - grid-scale only              
                
                !--in the boundary layer, add Wstar
                if(k >= kpbli(i,j) .and. k < LM)  wupdraft = wupdraft+zws(i,j) 

                IF(wupdraft > 0.1 .AND. wupdraft < 100.) THEN 

                ni   (1:n_modes)    =   max(AeroProps(i,j,k)%num(1:n_modes)*air_den,  zero_par)  ! unit: [m-3]
                rg   (1:n_modes)    =   max(AeroProps(i,j,k)%dpg(1:n_modes)*0.5*1.e6, zero_par)  ! unit: [um]
                sig0 (1:n_modes)    =       AeroProps(i,j,k)%sig(1:n_modes)                      ! unit: [um]
                bibar(1:n_modes)    =   max(AeroProps(i,j,k)%kap(1:n_modes),          zero_par)                 
              
                IF( tk >= 245.0) then   
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
                     
                     numbinit     = 0.
                     NACTL(i,j,k) = 0.
                     DO n=1,n_modes
                      numbinit     = numbinit    + AeroProps(i,j,k)%num(n)*air_den
                      NACTL(i,j,k) = NACTL(i,j,k)+ nact(n)
                     ENDDO
                     NACTL(i,j,k) = MIN(NACTL(i,j,k),0.99*numbinit)
                ENDIF ! tk>245
               ENDIF   ! updraft > 0.1
               ENDIF   ! plo > 34000.0
              

               IF( tk <= 268.0) then
                IF( (QC >= 0.5) .and. (QL >= 0.5)) then
                      ! Number of activated IN following deMott (2010) [#/m3]  
                         NACTI(i,j,k) = (1.e+3*ai*((273.16-tk)**bi) *  (naer_cb(i,j))**(ci*(273.16-tk)+di))  !#/m3
                ELSE
                      ! Number of activated IN following Wyser  
              
                 WC    = air_den*QC  !kg/m3
                 if (WC >= tiny(1.0)) then
                    BB = -2. + log10(1000.*WC/50.)*(1.e-3*(273.15-tk)**1.5)
                 else
                    BB = -6.
                 end if
                 BB    = MIN((MAX(BB,-6.)),-2.)  

                 RAUX  = 377.4 + 203.3 * BB+ 37.91 * BB **2 + 2.3696 * BB **3
                 RAUX  = (betai + (gamai + deltai * RAUX**3)**0.5)**0.33333
                 NACTI(i,j,k) = (3.* WC)/(4.*MAPL_PI*densic*(1.D-6*RAUX)**3)  !#/m3
               
                ENDIF  !Mixed phase
               ENDIF ! tk<=268
               !
               !
               !-- fix limit for NACTL/NACTI
               IF(NACTL(i,j,k) < NN_MIN) NACTL(i,j,k) = NN_MIN
               IF(NACTI(i,j,k) < NN_MIN) NACTI(i,j,k) = NN_MIN

        ENDDO;ENDDO;ENDDO

        deallocate(   rg, __STAT__)
        deallocate(   ni, __STAT__)
        deallocate(bibar, __STAT__)
        deallocate( nact, __STAT__)

      end if

      end if

      END SUBROUTINE Aer_Activation
      
!>----------------------------------------------------------------------------------------------------------------------
!!    12-12-06, DLW: Routine to set up the call to subr. ACTFRAC_MAT to calculate the 
!!                   activated fraction of the number and mass concentrations, 
!!                    as well as the number and mass concentrations activated 
!!                    for each of nmodes modes. The minimum dry radius for activation 
!!                    for each mode is also returned. 
!!
!!     Each mode is assumed to potentially contains 5 chemical species:
!!         (1) sulfate 
!!         (2) BC 
!!         (3) OC
!!         (4) mineral dust
!!         (5) sea salt 
!!
!!     The aerosol activation parameterizations are described in 
!!
!!         1. Abdul-Razzak et al.   1998, JGR, vol.103, p.6123-6131.
!!         2. Abdul-Razzak and Ghan 2000, JGR, vol.105, p.6837-6844. 
!!
!!     and values for many of the required parameters were taken from 
!!
!!         3. Ghan et al. 2001, JGR vol 106, p.5295-5316.
!!
!!     With the density of sea salt set to the value used in ref. 3 (1900 kg/m^3), this routine 
!!     yields values for the hygroscopicity parameters Bi in agreement with ref. 3. 
!!----------------------------------------------------------------------------------------------------------------------
      subroutine GetActFrac(nmodes  & !nmodes                                     &
      ,xnap                         & !ni              (1:nmodes)              & 
      ,rg                           & !0.5d+00*dgn_dry (1:nmodes)              & 
      ,sigmag                       & !sig0             (1:nmodes)              & 
      ,tkelvin                      & !tk              (i,j,k)                     &
      ,ptot                         & !pres             (i,j,k)                     & 
      ,wupdraft                     & !wupdraft             (i,j,k)                     & 
      ,nact                         & !nact             (i,j,k,1:nmodes)             &
      ,bibar)
      
      IMPLICIT NONE

      ! arguments.
      
      integer :: nmodes               !< number of modes [1]      
      real(AER_PR) :: xnap(nmodes)         !< number concentration for each mode [#/m^3]
      real(AER_PR) :: rg(nmodes)           !< geometric mean dry radius for each mode [um]
      real(AER_PR) :: sigmag(nmodes)       !< geometric standard deviation for each mode [um]
      real(AER_PR) :: tkelvin              !< absolute temperature [k]
      real(AER_PR) :: ptot                 !< ambient pressure [pa]
      real(AER_PR) :: wupdraft             !< updraft velocity [m/s]
!     real(AER_PR) :: ac(nmodes)           !< minimum dry radius for activation for each mode [um]
!     real(AER_PR) :: fracactn(nmodes)     !< activating fraction of number conc. for each mode [1]
      real(AER_PR) :: nact(nmodes)         !< activating number concentration for each mode [#/m^3]
      real(AER_PR) :: bibar(nmodes)        ! hygroscopicity parameter for each mode [1]

      ! local variables. 

      integer :: i, j                 ! loop counters       

      !--------------------------------------------------------------------------------------------------------------
      ! calculate the droplet activation parameters for each mode. 
      !--------------------------------------------------------------------------------------------------------------
      call ActFrac_Mat(nmodes,xnap,rg,sigmag,bibar,tkelvin,ptot,wupdraft,nact)

      end subroutine GetActFrac


!>----------------------------------------------------------------------------------------------------------------------
!!     12-12-06, DLW: Routine to calculate the activated fraction of the number 
!!                    and mass concentrations, as well as the number and mass 
!!                    concentrations activated for each of nmodes modes. The 
!!                    minimum dry radius for activation for each mode is also returned. 
!!
!!     The aerosol activation parameterizations are described in 
!!
!!         1. Abdul-Razzak et al.   1998, JGR, vol.103, p.6123-6131.
!!         2. Abdul-Razzak and Ghan 2000, JGR, vol.105, p.6837-6844. 
!! 
!!     This routine is for the multiple-aerosol type parameterization. 
!!----------------------------------------------------------------------------------------------------------------------
      subroutine ActFrac_Mat(nmodes,xnap,rg,sigmag,bibar,tkelvin,ptot,wupdraft,nact)

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
      real(AER_PR) :: ac(nmodes)        !< minimum dry radius for activation for each mode [um]
      real(AER_PR) :: fracactn(nmodes)  !< activating fraction of number conc. for each mode [1]
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
      real(AER_PR)            :: eta(nmodes)                    ! model parameter [1]  
      real(AER_PR)            :: zeta                           ! model parameter [1]  
      real(AER_PR)            :: xlogsigm(nmodes)               ! ln(sigmag) [1]   
      real(AER_PR)            :: a                              ! [m]
      real(AER_PR)            :: g                              ! [m^2/s]   
      real(AER_PR)            :: rdrp                           ! [m]   
      real(AER_PR)            :: f1                             ! [1]   
      real(AER_PR)            :: f2                             ! [1]
      real(AER_PR)            :: alpha                          ! [1/m]
      real(AER_PR)            :: gamma                          ! [m^3/kg]   
      real(AER_PR)            :: sm(nmodes)                     ! [1]   
      real(AER_PR)            :: dum                            ! [1/m]    
      real(AER_PR)            :: u                              ! argument to error function [1]
      real(AER_PR)            :: erf                            ! error function [1], but not declared in an f90 module 
      real(AER_PR)            :: smax                           ! maximum supersaturation [1]

      real(AER_PR)            :: aux1,aux2                           ! aux
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
      xlogsigm(:) = log(sigmag(:))                                                    ! [1] 
      smax = 0.0d+00                                                                  ! [1]
      
      do i=1, nmodes
    
        sm(i) = ( 2.0d+00/sqrt(bibar(i)) ) * ( a/(3.0*rg(i)) )**1.5d+00               ! [1] 
        eta(i) = dum**3 / (twopi*denh2o*gamma*xnap(i))                                ! [1] 

        !--------------------------------------------------------------------------------------------------------------
        ! write(1,'(a27,i4,4d15.5)')'i,eta(i),sm(i) =',i,eta(i),sm(i)
        !--------------------------------------------------------------------------------------------------------------
        f1 = 0.5d+00 * exp(2.50d+00 * xlogsigm(i)**2)                                 ! [1] 
        f2 = 1.0d+00 +     0.25d+00 * xlogsigm(i)                                     ! [1] 
        smax = smax + (   f1*(  zeta  / eta(i)              )**1.50d+00 &
                        + f2*(sm(i)**2/(eta(i)+3.0d+00*zeta))**0.75d+00 ) / sm(i)**2  ! [1] - eq. (6)
      enddo 
      smax = 1.0d+00 / sqrt(smax)                                                     ! [1]
      
      
      !aux1=0.0D0; aux2=0.0D0
      do i=1, nmodes

        ac(i)       = rg(i) * ( sm(i) / smax )**0.66666666666666667d+00               ! [um]

        u           = log(ac(i)/rg(i)) / ( sqrt2 * xlogsigm(i) )                      ! [1]
        fracactn(i) = 0.5d+00 * (1.0d+00 - erf(u))                                    ! [1]
        nact(i)     = fracactn(i) * xnap(i)                                           ! [#/m^3]
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


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine GcfMatrix(gammcf,a,x,gln)

      implicit none
      integer, parameter :: itmax=10000
      real(AER_PR), parameter :: eps=3.0d-07
      real(AER_PR), parameter :: fpmin=1.0d-30
      real(AER_PR) :: a,gammcf,gln,x
      integer :: i
      real(AER_PR) :: an,b,c,d,del,h
      !real(AER_PR) :: gammln   ! function names not declared in an f90 module 
      gln=gammln(a)
      b=x+1.0d+00-a
      c=1.0d+00/fpmin
      d=1.0d+00/b
      h=d
      do i=1,itmax
        an=-i*(i-a)
        b=b+2.0d+00
        d=an*d+b
        if(abs(d).lt.fpmin)d=fpmin
        c=b+an/c
        if(abs(c).lt.fpmin)c=fpmin
        d=1.0d+00/d
        del=d*c
        h=h*del
        if(abs(del-1.0d+00).lt.eps)goto 1
      enddo
      write(*,*)'AERO_ACTV: SUBROUTINE GCF: A TOO LARGE, ITMAX TOO SMALL', gammcf,a,x,gln
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      end subroutine GcfMatrix


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Gser(gamser,a,x,gln)

      implicit none
      integer, parameter :: itmax=10000  ! was itmax=100   in press et al. 
      real(AER_PR), parameter :: eps=3.0d-09  ! was eps=3.0d-07 in press et al.
      real(AER_PR) :: a,gamser,gln,x
      integer :: n
      real(AER_PR) :: ap,del,sum
      !real(AER_PR) :: gammln   ! function names not declared in an f90 module 
      gln=gammln(a)
      if(x.le.0.d+00)then
        if(x.lt.0.)stop 'aero_actv: subroutine gser: x < 0 in gser'
        gamser=0.d+00
        return
      endif
      ap=a
      sum=1.d+00/a
      del=sum
      do n=1,itmax
        ap=ap+1.d+00
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*eps)goto 1
      enddo
      write(*,*)'aero_actv: subroutine gser: a too large, itmax too small'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      end subroutine Gser


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      double precision function GammLn(xx)

      implicit none
      real(AER_PR) :: xx
      integer j
      double precision ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0,         &
      24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
      -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      enddo
      gammln=tmp+log(stp*ser/x)
      return
      end function GammLn 


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      double precision function Erf(x)
      implicit none
      real(AER_PR) :: x
!u    uses gammp
      !LFR  real(AER_PR) :: gammp   ! function names not declared in an f90 module 
      erf = 0.d0
      if(x.lt.0.0d+00)then
        erf=-gammp(0.5d0,x**2)
      else
        erf= gammp(0.5d0,x**2)
      endif
      return
      end function Erf


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      double precision function GammP(a,x)
      implicit none
      real(AER_PR) :: a,x
      real(AER_PR) :: gammcf,gamser,gln
      if(x.lt.0.0d+00.or.a.le.0.0d+00)then
        write(*,*)'aero_actv: function gammp: bad arguments'
      endif

      if(x.lt.a+1.0d+00)then
        call Gser(gamser,a,x,gln)
        gammp=gamser
      else
        call GcfMatrix(gammcf,a,x,gln)
        gammp=1.0d+00-gammcf
      endif
      return
      end function GammP
!>-----------------------------------------------------------------------------------------------------------------------


END MODULE Aer_Actv_Single_Moment

