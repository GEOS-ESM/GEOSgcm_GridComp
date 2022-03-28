MODULE Aer_Actv_Single_Moment
!
      USE aer_cloud, only: AerProps
      USE MAPL_ConstantsMod, only: MAPL_PI  
!-------------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      PUBLIC ::  Aer_Actv_1M_interface, USE_AEROSOL_NN
      PRIVATE
       real*8, parameter :: zero_par  =  1.0D-20
       real*8, parameter :: ai        =  0.0000594D0
       real*8, parameter :: bi        =  3.33D0
       real*8, parameter :: ci        =  0.0264D0
       real*8, parameter :: di        =  0.0033D0

       real*8, parameter :: betai     = -2.262D+3
       real*8, parameter :: gamai     =  5.113D+6
       real*8, parameter :: deltai    =  2.809D+3
       real*8, parameter :: densic    =  500.D0   !Ice crystal density in kgm-3

    
      !-- try 50.e6 over ocean, 300e6 over land   
       real, parameter :: NN_LAND     = 150.0e6
       real, parameter :: NN_OCEAN    =  30.0e6
       real, parameter :: NN_MIN      = NN_OCEAN

       LOGICAL  :: USE_AEROSOL_NN
      CONTAINS          

!>----------------------------------------------------------------------------------------------------------------------
!>----------------------------------------------------------------------------------------------------------------------

      SUBROUTINE Aer_Actv_1M_interface(IM,JM,LM,nmodes, t,plo,zlo,zle, qlcn, qicn, qlls, qils &
                                      , kpbl,zws, omega,FRLAND ,AeroProps, NACTL,NACTI)
      IMPLICIT NONE
      integer, intent(in)::IM,JM,LM,nmodes
      TYPE(AerProps), dimension (IM,JM,LM),intent(in )  :: AeroProps
      real, dimension (IM,JM,LM)  ,intent(in ) :: t,plo,omega,zlo, qlcn, qicn, qlls, qils
      real, dimension (IM,JM,0:LM),intent(in ) :: zle
      real, dimension (IM,JM)     ,intent(in ) :: zws,FRLAND
      integer, dimension (IM,JM)  ,intent(in ) :: kpbl     
 
      real, dimension (IM,JM,LM),intent(OUT) :: NACTL,NACTI
      
      real(8), dimension (IM,JM,LM,nmodes) ::nact 
      real(8), dimension (nmodes) :: sig0,rg,ni,bibar 
      real(8), dimension (IM,JM)  :: naer_cb
      real(8)                     :: wupdraft,tk,press,air_den,QC,QL,WC,BB,RAUX

      REAL :: numbinit
      integer :: i,j,k,n,kcb

      !--- activated aerosol # concentration for liq/ice phases (units: m^-3)
      numbinit = 0.
      NACTL    = 0.
      NACTI    = 0.
      WC       = 0.
      BB       = 0.
      RAUX     = 0.
            
      !--- determing aerosol number concentration at cloud base
      DO j=1,JM
        Do i=1,IM 
        !------check this
             kcb = kpbl(i,j)
        !------check this
             
             naer_cb(i,j) = zero_par
             !DO K=LM,kcb,-1
             k=min(kcb-1,LM-1); IF(K==0) stop "K==0- aer-act"
             tk           = DBLE(T(i,j,k))           ! K
             press        = DBLE(plo(i,j,k))*100.D0  ! Pa     
             air_den      = press*28.8d-3/8.31d0/tk  ! kg/m3
             DO n=1,nmodes
                if (AeroProps(i,j,k)%dpg(n) .ge. 0.5e-6) &
                naer_cb(i,j)= naer_cb(i,j) + DBLE(AeroProps(i,j,k)%num(n)) * air_den         
                naer_cb(i,j)= naer_cb(i,j)* 1.d+6  ! #/cm3
                naer_cb(i,j)= max(1.0d-1,min(naer_cb(i,j),100.0))
                                
             ENDDO
      ENDDO;ENDDO!;ENDDO
     
      DO k=LM,1,-1
       DO j=1,JM
        Do i=1,IM
              
              tk                 = DBLE(T(i,j,k))                          ! K
              press                 = DBLE(plo(i,j,k))*100.D0                  ! Pa   
              air_den                 = press*28.8d-3/8.31d0/tk                  ! kg/m3
              qc                 = DBLE(qicn(i,j,k)+qils(i,j,k))*1.D+3    ! g/kg
              ql                 = DBLE(qlcn(i,j,k)+qlls(i,j,k))*1.D+3    ! g/kg
              
              IF( plo(i,j,k)*100.0 > 34000.0) THEN 
                
                wupdraft           = -9.81*air_den*DBLE(omega(i,j,k))          ! m/s - grid-scale only              
                
                !--in the boundary layer, add Wstar
                if(k >= kpbl(i,j) .and. k < LM)  wupdraft = wupdraft+DBLE(zws(i,j)) 

                IF(wupdraft > 0.1 .AND. wupdraft < 100.) THEN 

                ni   (1:nmodes)    =   max(DBLE(AeroProps(i,j,k)%num(1:nmodes))*air_den,1.D-20)        ! unit: [m-3]
                rg   (1:nmodes)    =   1.D+6*max(1.0D-10,DBLE(AeroProps(i,j,k)%dpg(1:nmodes))*0.5d+00) ! unit: [um]
                sig0 (1:nmodes)    =   DBLE(AeroProps(i,j,k)%sig(1:nmodes))                            ! unit: [um]
                bibar(1:nmodes)    =   MAX(0.00001D0,DBLE(AeroProps(i,j,k)%kap(1:nmodes)))                 
              
                IF( tk >= 245.0D0) then   
                     call GetActFrac(nmodes                                   &
                                     ,ni              (1:nmodes)              &  
                                     ,rg              (1:nmodes)              & 
                                     ,sig0            (1:nmodes)              &  
                                     ,tk                                      &
                                     ,press                                   & 
                                     ,wupdraft                                & 
                                     ,nact            (i,j,k,1:nmodes)        &
                                     ,bibar           (1:nmodes)              &
                                     )
                     
                     numbinit     = 0.
                     NACTL(i,j,k) = 0.
                     DO n=1,nmodes
                      numbinit     = numbinit    + AeroProps(i,j,k)%num(n)*air_den
                      NACTL(i,j,k) = NACTL(i,j,k)+ real(nact(i,j,k,n),4) 
                     ENDDO
                
                 
                     NACTL(i,j,k) = MIN(NACTL(i,j,k),0.99*numbinit)
                
 
                ENDIF ! tk>245
               ENDIF   ! updraft > 0.1
               ENDIF   ! 100*plo > 34000.0
              

               IF( tk <= 268.0D0) then
                IF( (QC >= 0.5) .and. (QL >= 0.5)) then
                      ! Number of activated IN following deMott (2010) [#/m3]  
                         NACTI(i,j,k) = real((1.e+3*ai*((273.16-tk)**bi) *  (naer_cb(i,j))**(ci*(273.16-tk)+di)),4)  !#/m3
                  ELSE   !tk<243
                      ! Number of activated IN following Wyser  
              
                 WC    = air_den*QC  !kg/m3
                 if (WC >= tiny(1.0d0)) then
                    BB     =  -2. + log10(1000.*WC/50.)*(1.e-3*(273.15-tk)**1.5)
                 else
                    BB = -6.
                 end if
                 BB     = MIN((MAX(BB,-6.)),-2.)  

                 RAUX   = 377.4 + 203.3 * BB+ 37.91 * BB **2 + 2.3696 * BB **3
                 RAUX   = (betai + (gamai + deltai * RAUX**3)**0.5)**0.33333
                 NACTI(i,j,k) = real((3.* WC)/(4.*MAPL_PI*densic*(1.D-6*RAUX)**3),4)  !#/m3
               
                ENDIF  !Mixed phase

               ENDIF ! tk<=268
               !
               !
               !-- fix limit for NACTL/NACTI
               IF(NACTL(i,j,k) < NN_MIN) NACTL(i,j,k) = FRLAND(i,j)*NN_LAND + (1.0-FRLAND(i,j))*NN_OCEAN

        ENDDO;ENDDO;ENDDO

      END SUBROUTINE Aer_Actv_1M_interface
      
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
      real(8) :: xnap(nmodes)         !< number concentration for each mode [#/m^3]
      real(8) :: rg(nmodes)           !< geometric mean dry radius for each mode [um]
      real(8) :: sigmag(nmodes)       !< geometric standard deviation for each mode [um]
      real(8) :: tkelvin              !< absolute temperature [k]
      real(8) :: ptot                 !< ambient pressure [pa]
      real(8) :: wupdraft             !< updraft velocity [m/s]
!     real(8) :: ac(nmodes)           !< minimum dry radius for activation for each mode [um]
!     real(8) :: fracactn(nmodes)     !< activating fraction of number conc. for each mode [1]
      real(8) :: nact(nmodes)         !< activating number concentration for each mode [#/m^3]
      real(8) :: bibar(nmodes)        ! hygroscopicity parameter for each mode [1]

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
      real(8) :: xnap(nmodes)      !< number concentration for each mode [#/m^3]
!     real(8) :: xmap(nmodes)      !< mass   concentration for each mode [ug/m^3]
      real(8) :: rg(nmodes)        !< geometric mean radius for each mode [um]
      real(8) :: sigmag(nmodes)    !< geometric standard deviation for each mode [um]
      real(8) :: bibar(nmodes)     !< hygroscopicity parameter for each mode [1]
      real(8) :: tkelvin           !< absolute temperature [k]
      real(8) :: ptot              !< ambient pressure [pa]
      real(8) :: wupdraft          !< updraft velocity [m/s]
      real(8) :: ac(nmodes)        !< minimum dry radius for activation for each mode [um]
      real(8) :: fracactn(nmodes)  !< activating fraction of number conc. for each mode [1]
      real(8) :: nact(nmodes)      !< activating number concentration for each mode [#/m^3]

      ! parameters.
      
      real(8), parameter :: pi            = 3.141592653589793d+00
      real(8), parameter :: twopi         = 2.0d+00 * pi
      real(8), parameter :: sqrt2         = 1.414213562d+00
      real(8), parameter :: threesqrt2by2 = 1.5d+00 * sqrt2

      real(8), parameter :: avgnum   = 6.0221367d+23       ! [1/mol]
      real(8), parameter :: rgasjmol = 8.31451d+00         ! [j/mol/k]
      real(8), parameter :: wmolmass = 18.01528d-03        ! molar mass of h2o     [kg/mol]
      real(8), parameter :: amolmass = 28.966d-03          ! molar mass of air     [kg/mol]
      real(8), parameter :: asmolmss = 132.1406d-03        ! molar mass of nh42so4 [kg/mol]
      real(8), parameter :: denh2o   = 1.00d+03            ! density of water [kg/m^3]
      real(8), parameter :: denamsul = 1.77d+03            ! density of pure ammonium sulfate [kg/m^3]
      real(8), parameter :: xnuamsul = 3.00d+00            ! # of ions formed when the salt is dissolved in water [1]
      real(8), parameter :: phiamsul = 1.000d+00           ! osmotic coefficient value in a-r 1998. [1] 
      real(8), parameter :: gravity  = 9.81d+00            ! grav. accel. at the earth's surface [m/s/s] 
      real(8), parameter :: heatvap  = 40.66d+03/wmolmass  ! latent heat of vap. for water and tnbp [j/kg] 
      real(8), parameter :: cpair    = 1006.0d+00          ! heat capacity of air [j/kg/k] 
      real(8), parameter :: t0dij    = 273.15d+00          ! reference temp. for dv [k] 
      real(8), parameter :: p0dij    = 101325.0d+00        ! reference pressure for dv [pa] 
      real(8), parameter :: dijh2o0  = 0.211d-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
      !----------------------------------------------------------------------------------------------------------------    
      ! real(8), parameter :: t0dij    = 283.15d+00          ! reference temp. for dv [k] 
      ! real(8), parameter :: p0dij    = 80000.0d+00         ! reference pressure for dv [pa] 
      ! real(8), parameter :: dijh2o0  = 0.300d-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
      !----------------------------------------------------------------------------------------------------------------
      real(8), parameter :: deltav   = 1.096d-07           ! vapor jump length [m]  
      real(8), parameter :: deltat   = 2.160d-07           ! thermal jump length [m]  
      real(8), parameter :: alphac   = 1.000d+00           ! condensation mass accommodation coefficient [1]  
      real(8), parameter :: alphat   = 0.960d+00           ! thermal accommodation coefficient [1]  

      ! local variables. 

      integer            :: i                              ! loop counter 
      real(8)            :: dv                             ! diffusion coefficient for water [m^2/s] 
      real(8)            :: dvprime                        ! modified diffusion coefficient for water [m^2/s] 
      real(8)            :: dumw, duma                     ! scratch variables [s/m] 
      real(8)            :: wpe                            ! saturation vapor pressure of water [pa]  
      real(8)            :: surten                         ! surface tension of air-water interface [j/m^2] 
      real(8)            :: xka                            ! thermal conductivity of air [j/m/s/k]  
      real(8)            :: xkaprime                       ! modified thermal conductivity of air [j/m/s/k]  
      real(8)            :: eta(nmodes)                    ! model parameter [1]  
      real(8)            :: zeta                           ! model parameter [1]  
      real(8)            :: xlogsigm(nmodes)               ! ln(sigmag) [1]   
      real(8)            :: a                              ! [m]
      real(8)            :: g                              ! [m^2/s]   
      real(8)            :: rdrp                           ! [m]   
      real(8)            :: f1                             ! [1]   
      real(8)            :: f2                             ! [1]
      real(8)            :: alpha                          ! [1/m]
      real(8)            :: gamma                          ! [m^3/kg]   
      real(8)            :: sm(nmodes)                     ! [1]   
      real(8)            :: dum                            ! [1/m]    
      real(8)            :: u                              ! argument to error function [1]
      real(8)            :: erf                            ! error function [1], but not declared in an f90 module 
      real(8)            :: smax                           ! maximum supersaturation [1]

      real(8)            :: aux1,aux2                           ! aux
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
      real(8), parameter :: eps=3.0d-07
      real(8), parameter :: fpmin=1.0d-30
      real(8) :: a,gammcf,gln,x
      integer :: i
      real(8) :: an,b,c,d,del,h
      !real(8) :: gammln   ! function names not declared in an f90 module 
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
      real(8), parameter :: eps=3.0d-09  ! was eps=3.0d-07 in press et al.
      real(8) :: a,gamser,gln,x
      integer :: n
      real(8) :: ap,del,sum
      !real(8) :: gammln   ! function names not declared in an f90 module 
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
      real(8) :: xx
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
      real(8) :: x
!u    uses gammp
      !LFR  real(8) :: gammp   ! function names not declared in an f90 module 
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
      real(8) :: a,x
      real(8) :: gammcf,gamser,gln
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

