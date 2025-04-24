MODULE Aer_Actv_Single_Moment
!
#include "MAPL_Generic.h"

      USE ESMF
      USE MAPL
      USE aer_cloud, only: AerPropsNew
!-------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      PUBLIC ::  Aer_Activation, USE_BERGERON, USE_AEROSOL_NN, R_AIR
      PRIVATE

       ! Real kind for activation.
       integer,public,parameter :: AER_PR = MAPL_R4

       real        , parameter :: R_AIR     =  3.47e-3 !m3 Pa kg-1K-1
       real(AER_PR), parameter :: zero_par  =  tiny(1.0)   ! small non-zero value
       real(AER_PR), parameter :: ai        =  0.0000594
       real(AER_PR), parameter :: bi        =  3.33
       real(AER_PR), parameter :: ci        =  0.0264
       real(AER_PR), parameter :: di        =  0.0033

       real(AER_PR), parameter :: betai     = -2.262e+3
       real(AER_PR), parameter :: gamai     =  5.113e+6
       real(AER_PR), parameter :: deltai    =  2.809e+3
       real(AER_PR), parameter :: densic    =  917.0   !Ice crystal density in kgm-3

       real, parameter :: NN_MIN      =  100.0e6
       real, parameter :: NN_MAX      =  500.0e6

       LOGICAL  :: USE_BERGERON = .TRUE.
       LOGICAL  :: USE_AEROSOL_NN = .TRUE.
      CONTAINS          

!>----------------------------------------------------------------------------------------------------------------------
!>----------------------------------------------------------------------------------------------------------------------

      SUBROUTINE Aer_Activation(MAPL, IM,JM,LM, q, t, plo, ple, tke, vvel, FRLAND, &
                                AeroPropsNew, aero_aci, NACTL, NACTI, NWFA,  &
                                NN_LAND, NN_OCEAN, need_extra_fields)
      IMPLICIT NONE
      type (MAPL_MetaComp), pointer   :: MAPL
      integer, intent(in)::IM,JM,LM
      TYPE(AerPropsNew), dimension (:), intent(inout) :: AeroPropsNew
      type(ESMF_State)            ,intent(inout) :: aero_aci
      real, dimension (IM,JM,LM)  ,intent(in ) :: plo ! Pa
      real, dimension (IM,JM,0:LM),intent(in ) :: ple ! Pa
      real, dimension (IM,JM,LM)  ,intent(in ) :: q,t,tke,vvel
      real, dimension (IM,JM)     ,intent(in ) :: FRLAND
      real                        ,intent(in ) :: NN_LAND, NN_OCEAN     
      logical                     ,intent(in ) :: need_extra_fields 
 
      real, dimension (IM,JM,LM),intent(OUT) :: NACTL, NACTI, NWFA
      
      real(AER_PR), allocatable, dimension (:,:,:) :: sig0,rg,ni,bibar,nact 
      real(AER_PR), dimension(IM,JM)               :: wupdraft,tk,press,air_den

      integer, parameter :: ALT_MAXSTR=64
      character(len=ALT_MAXSTR)      :: aci_field_name
      real, pointer, dimension(:,:)   :: aci_ptr_2d
      real, pointer, dimension(:,:,:) :: aci_ptr_3d
      character(len=ALT_MAXSTR), allocatable, dimension(:) :: aero_aci_modes
      integer                         :: ACI_STATUS

      integer :: n_modes
      REAL :: numbinit(IM,JM)
      integer :: k,n,rc
      integer :: nn

      character(len=ESMF_MAXSTR)              :: IAm="Aer_Activation"
      integer                                 :: STATUS

      NWFA = 0.0

      if (.not. USE_AEROSOL_NN) then

        do k = 1, LM
          NACTL(:,:,k) = NN_LAND*FRLAND + NN_OCEAN*(1.0-FRLAND)
          NACTI(:,:,k) = NN_LAND*FRLAND + NN_OCEAN*(1.0-FRLAND)
        end do

        return
     end if

     call ESMF_AttributeGet(aero_aci, name='number_of_aerosol_modes', value=n_modes, __RC__)

     if (n_modes == 0) return

     call MAPL_TimerOn (MAPL,"----AERO_ACTIVATE_1",__RC__)
     
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

     ACTIVATION_PROPERTIES: do n = 1, n_modes
        call ESMF_AttributeSet(aero_aci, name='aerosol_mode', value=trim(aero_aci_modes(n)), __RC__)
        ! call WRITE_PARALLEL (trim(aero_aci_modes(n)))  

        ! execute the aerosol activation properties method 
        call ESMF_MethodExecute(aero_aci, label='aerosol_activation_properties', userRC=ACI_STATUS, RC=STATUS)
        VERIFY_(ACI_STATUS)
        VERIFY_(STATUS)

        ! copy out aerosol activation properties
        call ESMF_AttributeGet(aero_aci, name='aerosol_number_concentration', value=aci_field_name, __RC__)
        call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
        AeroPropsNew(n)%num = aci_ptr_3d

        call ESMF_AttributeGet(aero_aci, name='aerosol_dry_size', value=aci_field_name, __RC__)
        call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
        AeroPropsNew(n)%dpg = aci_ptr_3d
        ! if (MAPL_am_I_root()) print *, AeroPropsNew(n)%dpg(1,1,1)

        call ESMF_AttributeGet(aero_aci, name='width_of_aerosol_mode', value=aci_field_name, __RC__)
        call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
        AeroPropsNew(n)%sig = aci_ptr_3d

        call ESMF_AttributeGet(aero_aci, name='aerosol_hygroscopicity', value=aci_field_name, __RC__)
        call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
        AeroPropsNew(n)%kap = aci_ptr_3d
        ! if (MAPL_am_I_root()) print *, AeroPropsNew(n)%kap(1,1,1)

        if (need_extra_fields) then

           call ESMF_AttributeGet(aero_aci, name='aerosol_density', value=aci_field_name, __RC__)
           call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
           AeroPropsNew(n)%den = aci_ptr_3d

           call ESMF_AttributeGet(aero_aci, name='fraction_of_dust_aerosol', value=aci_field_name, __RC__)
           call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
           AeroPropsNew(n)%fdust = aci_ptr_3d

           call ESMF_AttributeGet(aero_aci, name='fraction_of_soot_aerosol', value=aci_field_name, __RC__)
           call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
           AeroPropsNew(n)%fsoot = aci_ptr_3d

           call ESMF_AttributeGet(aero_aci, name='fraction_of_organic_aerosol', value=aci_field_name, __RC__)
           call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
           AeroPropsNew(n)%forg = aci_ptr_3d

        endif

        AeroPropsNew(n)%nmods = n_modes

        where (AeroPropsNew(n)%kap > 0.4)
           NWFA = NWFA + AeroPropsNew(n)%num
        end where

     end do ACTIVATION_PROPERTIES

     deallocate(aero_aci_modes, __STAT__)

     call MAPL_Timeroff(MAPL,"----AERO_ACTIVATE_1",__RC__)

     call MAPL_TimerOn (MAPL,"----AERO_ACTIVATE_2",__RC__)
     !--- activated aerosol # concentration for liq/ice phases (units: m^-3)

     allocate( sig0(IM,JM,n_modes), __STAT__)
     allocate(   rg(IM,JM,n_modes), __STAT__)
     allocate(   ni(IM,JM,n_modes), __STAT__)
     allocate(bibar(IM,JM,n_modes), __STAT__)
     allocate( nact(IM,JM,n_modes), __STAT__)

!$OMP parallel do default(none) shared(IM,JM,LM,n_modes,T,plo,vvel,tke,MAPL_RGAS,zero_par, &
!$OMP                                  AeroPropsNew,NACTL,NACTI,NN_MIN,NN_MAX,ai,bi,ci,di) &
!$OMP                           private(k,n,tk,press,air_den,wupdraft,ni,rg,bibar,sig0,nact)
     DO k=1,LM

           tk                 = T(:,:,k)                         ! K
           press              = plo(:,:,k)                       ! Pa   
           air_den            = press/(MAPL_RGAS*tk)             ! kg/m3
           wupdraft           = max(zero_par,vvel(:,:,k) + SQRT(tke(:,:,k)))

           ! Liquid Clouds
           ni = tiny(1.0)
           DO n=1,n_modes
              where (AeroPropsNew(n)%kap(:,:,k) > 0.4) &
              ni   (:,:,n) =   max(AeroPropsNew(n)%num(:,:,k)*air_den,  zero_par)  ! unit: [m-3]
              rg   (:,:,n) =   max(AeroPropsNew(n)%dpg(:,:,k)*0.5e6,    zero_par)  ! unit: [um]
              bibar(:,:,n) =   max(AeroPropsNew(n)%kap(:,:,k),          zero_par)                 
              sig0 (:,:,n) =       AeroPropsNew(n)%sig(:,:,k)
           ENDDO
           call GetActFrac(IM*JM, n_modes    &
                ,      ni(:,:,1)   &
                ,      rg(:,:,1)   &
                ,    sig0(:,:,1)   &
                ,   bibar(:,:,1)   &
                ,      tk(:,:)     &
                ,   press(:,:)     &
                ,wupdraft(:,:)     &
                ,    nact(:,:,1)   &
                )
           numbinit = 0.
           NACTL(:,:,k) = 0.
           DO n=1,n_modes
              where (AeroPropsNew(n)%kap(:,:,k) > 0.4)
                 numbinit = numbinit + AeroPropsNew(n)%num(:,:,k)
                 NACTL(:,:,k)= NACTL(:,:,k) + nact(:,:,n) !#/m3
              end where
           ENDDO
           numbinit = numbinit * air_den ! #/m3
           NACTL(:,:,k) = MIN(NACTL(:,:,k),0.99*numbinit)
           NACTL(:,:,k) = MAX(MIN(NACTL(:,:,k),NN_MAX),NN_MIN)

           ! Ice Clouds
           numbinit = 0.
           DO n=1,n_modes
              where ( (AeroPropsNew(n)%dpg(:,:,k) .ge. 0.5e-6) .and. & ! diameters > 0.5 microns
                   (AeroPropsNew(n)%kap(:,:,k) .gt. 0.4) )
                 numbinit = numbinit + AeroPropsNew(n)%num(:,:,k)
              end where
           ENDDO
           numbinit = numbinit * air_den ! #/m3
           ! Number of activated IN following deMott (2010) [#/m3]
           NACTI(:,:,k) = (ai*(max(0.0,(MAPL_TICE-tk))**bi)) * (numbinit**(ci*max((MAPL_TICE-tk),0.0)+di)) !#/m3
           NACTI(:,:,k) = MAX(MIN(NACTI(:,:,k),NN_MAX),NN_MIN)

        ENDDO

        deallocate( sig0, __STAT__)
        deallocate(   rg, __STAT__)
        deallocate(   ni, __STAT__)
        deallocate(bibar, __STAT__)
        deallocate( nact, __STAT__)

     call MAPL_TimerOff(MAPL,"----AERO_ACTIVATE_2",__RC__)
      

      END SUBROUTINE Aer_Activation
      
!>----------------------------------------------------------------------------------------------------------------------
!!     12-12-06, DLW: Routine to calculate the activated fraction of the number 
!!                    and mass concentrations, as well as the number and mass 
!!                    concentrations activated for each of nmodes modes. The 
!!                    minimum dry radius for activation for each mode is also returned. 
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
!!
!!     The aerosol activation parameterizations are described in 
!!
!!         1. Abdul-Razzak et al.   1998, JGR, vol.103, p.6123-6131.
!!         2. Abdul-Razzak and Ghan 2000, JGR, vol.105, p.6837-6844. 
!! 
!!     This routine is for the multiple-aerosol type parameterization. 
!!----------------------------------------------------------------------------------------------------------------------
      subroutine GetActFrac(im, nmodes,xnap,rg,sigmag,bibar,tkelvin,ptot,wupdraft,nact)

      IMPLICIT NONE

      ! Arguments.

      integer :: im
      integer :: nmodes            !< number of modes [1]      
      real(AER_PR) :: xnap(im,nmodes)      !< number concentration for each mode [#/m^3]
!     real(AER_PR) :: xmap(im,nmodes)      !< mass   concentration for each mode [ug/m^3]
      real(AER_PR) :: rg(im,nmodes)        !< geometric mean radius for each mode [um]
      real(AER_PR) :: sigmag(im,nmodes)    !< geometric standard deviation for each mode [um]
      real(AER_PR) :: bibar(im,nmodes)     !< hygroscopicity parameter for each mode [1]
      real(AER_PR) :: tkelvin(im)           !< absolute temperature [k]
      real(AER_PR) :: ptot(im)              !< ambient pressure [pa]
      real(AER_PR) :: wupdraft(im)          !< updraft velocity [m/s]
      real(AER_PR) :: ac(im,nmodes)        !< minimum dry radius for activation for each mode [um]
      real(AER_PR) :: fracactn(im,nmodes)  !< activating fraction of number conc. for each mode [1]
      real(AER_PR) :: nact(im,nmodes)      !< activating number concentration for each mode [#/m^3]

      ! parameters.
      
      real(AER_PR), parameter :: pi            = 3.141592653589793
      real(AER_PR), parameter :: twopi         = 2.0 * pi
      real(AER_PR), parameter :: sqrt2         = 1.414213562
      real(AER_PR), parameter :: threesqrt2by2 = 1.5 * sqrt2

      real(AER_PR), parameter :: avgnum   = 6.0221367d+23       ! [1/mol]
      real(AER_PR), parameter :: rgasjmol = 8.31451         ! [j/mol/k]
      real(AER_PR), parameter :: wmolmass = 18.01528e-03        ! molar mass of h2o     [kg/mol]
      real(AER_PR), parameter :: amolmass = 28.966e-03          ! molar mass of air     [kg/mol]
      real(AER_PR), parameter :: asmolmss = 132.1406e-03        ! molar mass of nh42so4 [kg/mol]
      real(AER_PR), parameter :: denh2o   = 1.00d+03            ! density of water [kg/m^3]
      real(AER_PR), parameter :: denamsul = 1.77d+03            ! density of pure ammonium sulfate [kg/m^3]
      real(AER_PR), parameter :: xnuamsul = 3.00            ! # of ions formed when the salt is dissolved in water [1]
      real(AER_PR), parameter :: phiamsul = 1.000           ! osmotic coefficient value in a-r 1998. [1] 
      real(AER_PR), parameter :: gravity  = 9.81            ! grav. accel. at the earth's surface [m/s/s] 
      real(AER_PR), parameter :: heatvap  = 40.66d+03/wmolmass  ! latent heat of vap. for water and tnbp [j/kg] 
      real(AER_PR), parameter :: cpair    = 1006.0          ! heat capacity of air [j/kg/k] 
      real(AER_PR), parameter :: t0dij    = 273.15          ! reference temp. for dv [k] 
      real(AER_PR), parameter :: p0dij    = 101325.0        ! reference pressure for dv [pa] 
      real(AER_PR), parameter :: dijh2o0  = 0.211e-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
      !----------------------------------------------------------------------------------------------------------------    
      ! real(AER_PR), parameter :: t0dij    = 283.15          ! reference temp. for dv [k] 
      ! real(AER_PR), parameter :: p0dij    = 80000.0         ! reference pressure for dv [pa] 
      ! real(AER_PR), parameter :: dijh2o0  = 0.300e-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
      !----------------------------------------------------------------------------------------------------------------
      real(AER_PR), parameter :: deltav   = 1.096e-07           ! vapor jump length [m]  
      real(AER_PR), parameter :: deltat   = 2.160e-07           ! thermal jump length [m]  
      real(AER_PR), parameter :: alphac   = 1.000           ! condensation mass accommodation coefficient [1]  
      real(AER_PR), parameter :: alphat   = 0.960           ! thermal accommodation coefficient [1]  

      ! local variables. 

      integer            :: i, n                              ! loop counter 
      real(AER_PR)            :: dv(im)                             ! diffusion coefficient for water [m^2/s] 
      real(AER_PR)            :: dvprime(im)                        ! modified diffusion coefficient for water [m^2/s] 
      real(AER_PR)            :: dumw(im), duma(im)                     ! scratch variables [s/m] 
      real(AER_PR)            :: wpe(im)                            ! saturation vapor pressure of water [pa]  
      real(AER_PR)            :: surten(im)                         ! surface tension of air-water interface [j/m^2] 
      real(AER_PR)            :: xka(im)                            ! thermal conductivity of air [j/m/s/k]  
      real(AER_PR)            :: xkaprime(im)                       ! modified thermal conductivity of air [j/m/s/k]  
      real(AER_PR)            :: eta(im,nmodes)                    ! model parameter [1]  
      real(AER_PR)            :: zeta(im)                           ! model parameter [1]  
      real(AER_PR)            :: xlogsigm(im,nmodes)               ! ln(sigmag) [1]   
      real(AER_PR)            :: a(im)                              ! [m]
      real(AER_PR)            :: g(im)                              ! [m^2/s]   
      real(AER_PR)            :: rdrp(im)                           ! [m]   
      real(AER_PR)            :: f1(im)                             ! [1]   
      real(AER_PR)            :: f2(im)                             ! [1]
      real(AER_PR)            :: alpha(im)                          ! [1/m]
      real(AER_PR)            :: gamma(im)                          ! [m^3/kg]   
      real(AER_PR)            :: sm(im,nmodes)                     ! [1]   
      real(AER_PR)            :: dum(im)                            ! [1/m]    
      real(AER_PR)            :: u(im)                              ! argument to error function [1]
      real(AER_PR)            :: erf                            ! error function [1], but not declared in an f90 module 
      real(AER_PR)            :: smax(im)                           ! maximum supersaturation [1]

!----------------------------------------------------------------------------------------------------------------------
!     rdrp is the radius value used in eqs.(17) & (18) and was adjusted to yield eta and zeta 
!     values close to those given in a-z et al. 1998 figure 5. 
!----------------------------------------------------------------------------------------------------------------------
      rdrp = 0.105e-06   ! [m] tuned to approximate the results in figures 1-5 in a-z et al. 1998.  
!----------------------------------------------------------------------------------------------------------------------
!     these variables are common to all modes and need only be computed once. 
!----------------------------------------------------------------------------------------------------------------------
      dv = dijh2o0*(p0dij/ptot)*(tkelvin/t0dij)**1.94                 ! [m^2/s] (p&k,2nd ed., p.503)
      surten = 76.10e-03 - 0.155e-03 * (tkelvin-273.15)               ! [j/m^2] 
      wpe = exp( 77.34491296 - 7235.424651/tkelvin - 8.2*log(tkelvin) + tkelvin*5.7113e-03 )  ! [pa] 
      dumw = sqrt(twopi*wmolmass/rgasjmol/tkelvin)                        ! [s/m] 
      dvprime = dv / ( (rdrp/(rdrp+deltav)) + (dv*dumw/(rdrp*alphac)) )   ! [m^2/s] - eq. (17) 
      xka = (5.69+0.017*(tkelvin-273.15))*418.4e-05           ! [j/m/s/k] (0.0238 j/m/s/k at 273.15 k)
      duma = sqrt(twopi*amolmass/rgasjmol/tkelvin)                        ! [s/m]
      xkaprime = xka / ( ( rdrp/(rdrp+deltat) ) + ( xka*duma/(rdrp*alphat*denh2o*cpair) ) )   ! [j/m/s/k]
      g = 1.0 / ( (denh2o*rgasjmol*tkelvin) / (wpe*dvprime*wmolmass) &
                      + ( (heatvap*denh2o) / (xkaprime*tkelvin) ) &
                      * ( (heatvap*wmolmass) / (rgasjmol*tkelvin) - 1.0 ) )               ! [m^2/s]
      a = (2.0*surten*wmolmass)/(denh2o*rgasjmol*tkelvin)                                 ! [m] 
      alpha = (gravity/(rgasjmol*tkelvin))*((wmolmass*heatvap)/(cpair*tkelvin) - amolmass)    ! [1/m] 
      gamma = (rgasjmol*tkelvin)/(wpe*wmolmass) &
            + (wmolmass*heatvap*heatvap)/(cpair*ptot*amolmass*tkelvin)                        ! [m^3/kg]
      dum = sqrt(alpha*wupdraft/g)                  ! [1/m] 
      zeta = 2.*a*dum/3.                    ! [1] 
!----------------------------------------------------------------------------------------------------------------
      ! write(1,'(a27,4d15.5)')'surten,wpe,a            =',surten,wpe,a
      ! write(1,'(a27,4d15.5)')'xka,xkaprime,dv,dvprime =',xka,xkaprime,dv,dvprime
      ! write(1,'(a27,4d15.5)')'alpha,gamma,g, zeta     =',alpha,gamma,g,zeta
!----------------------------------------------------------------------------------------------------------------------
!     these variables must be computed for each mode. 
!----------------------------------------------------------------------------------------------------------------------
      xlogsigm(:,:) = log(sigmag(:,:))                                                    ! [1] 
      smax = 0.0                                                                  ! [1]

      do n=1, nmodes
    
        sm(:,n) = ( 2.0/sqrt(bibar(i,n)) ) * ( a/(3.0*rg(:,n)) )**1.5               ! [1] 
        eta(:,n) = dum**3 / (twopi*denh2o*gamma*xnap(:,n))                                ! [1] 

        !--------------------------------------------------------------------------------------------------------------
        ! write(1,'(a27,i4,4d15.5)')'i,eta(i),sm(i) =',i,eta(i),sm(i)
        !--------------------------------------------------------------------------------------------------------------
        f1 = 0.5 * exp(2.50 * xlogsigm(:,n)**2)                                 ! [1] 
        f2 = 1.0 +     0.25 * xlogsigm(:,n)                                     ! [1] 
        smax = smax + (   f1*(  zeta  / eta(:,n)              )**1.50 &
                        + f2*(sm(i,n)**2/(eta(:,n)+3.0*zeta))**0.75 ) / sm(:,n)**2  ! [1] - eq. (6)
     enddo
      smax = 1.0 / sqrt(smax)                                                     ! [1]

      do n=1, nmodes

        ac(:,n)       = rg(:,n) * ( sm(:,n) / smax )**0.66666666666666667               ! [um]

        u           = log(ac(:,n)/rg(:,n)) / ( sqrt2 * xlogsigm(:,n) )                      ! [1]
        fracactn(:,n) = 0.5 * (1.0 - erf(u))                                    ! [1]
        nact(:,n)     = min(fracactn(:,n),0.99) * xnap(:,n)                             ! [#/m^3]
        
        !if(fracactn(i) .gt. 0.9999999 ) then
        !   write(*,*)i,ac(i),u,fracactn(i),xnap(i)
        !   print*,' xxx',i,ac(i),u,fracactn(i),xnap(i)
        !   stop
        !endif

     end do

      return
      end subroutine GetActFrac


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine GcfMatrix(gammcf,a,x,gln)

      implicit none
      integer, parameter :: itmax=10000
      real(AER_PR), parameter :: eps=3.0e-07
      real(AER_PR), parameter :: fpmin=1.0e-30
      real(AER_PR) :: a,gammcf,gln,x
      integer :: i
      real(AER_PR) :: an,b,c,d,del,h
      gln=gammln(a)
      b=x+1.0-a
      c=1.0/fpmin
      d=1.0/b
      h=d
      do i=1,itmax
        an=-i*(i-a)
        b=b+2.0
        d=an*d+b
        if(abs(d).lt.fpmin)d=fpmin
        c=b+an/c
        if(abs(c).lt.fpmin)c=fpmin
        d=1.0/d
        del=d*c
        h=h*del
        if(abs(del-1.0).lt.eps)goto 1
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
      real(AER_PR), parameter :: eps=3.0e-09  ! was eps=3.0e-07 in press et al.
      real(AER_PR) :: a,gamser,gln,x
      integer :: n
      real(AER_PR) :: ap,del,sum
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)stop 'aero_actv: subroutine gser: x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do n=1,itmax
        ap=ap+1.
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
      real(AER_PR) function GammLn(xx)

      implicit none
      real(AER_PR) :: xx
      integer j
      real(AER_PR) ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146,-86.50532032941677,         &
      24.01409824083091,-1.231739572450155,.1208650973866179e-2, &
      -.5395239384953e-5,2.5066282746310005/
      x=xx
      y=x
      tmp=x+5.5
      tmp=(x+0.5)*log(tmp)-tmp
      ser=1.000000000190015
      do j=1,6
        y=y+1.
        ser=ser+cof(j)/y
      enddo
      gammln=tmp+log(stp*ser/x)
      return
      end function GammLn 


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      real(AER_PR) function Erf(x)
      implicit none
      real(AER_PR) :: x
      erf = 0.
      if(x.lt.0.0)then
        erf=-gammp(0.5,x**2)
      else
        erf= gammp(0.5,x**2)
      endif
      return
      end function Erf


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      real(AER_PR) function GammP(a,x)
      implicit none
      real(AER_PR) :: a,x
      real(AER_PR) :: gammcf,gamser,gln
      if(x.lt.0.0.or.a.le.0.0)then
        write(*,*)'aero_actv: function gammp: bad arguments'
      endif

      if(x.lt.a+1.0)then
        call Gser(gamser,a,x,gln)
        gammp=gamser
      else
        call GcfMatrix(gammcf,a,x,gln)
        gammp=1.0-gammcf
      endif
      return
      end function GammP
!>-----------------------------------------------------------------------------------------------------------------------


END MODULE Aer_Actv_Single_Moment

