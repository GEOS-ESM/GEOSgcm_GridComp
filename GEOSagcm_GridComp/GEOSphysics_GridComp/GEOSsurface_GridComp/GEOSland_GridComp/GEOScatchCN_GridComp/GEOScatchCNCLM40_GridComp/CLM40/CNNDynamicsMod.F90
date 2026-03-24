module CNNDynamicsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNNDynamicsMod
!
! !DESCRIPTION:
! Module for mineral nitrogen dynamics (deposition, fixation, leaching)
! for coupled carbon-nitrogen code.
!
! !USES:
   implicit none
   save
   private
! !PUBLIC MEMBER FUNCTIONS:
   public :: CNNDeposition
   public :: CNNFixation
   public :: CNNLeaching
!
! !REVISION HISTORY:
! 6/1/04: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNNDeposition
!
! !INTERFACE:
subroutine CNNDeposition( lbc, ubc )
!
! !DESCRIPTION:
! On the radiation time step, update the nitrogen deposition rate
! from atmospheric forcing. For now it is assumed that all the atmospheric
! N deposition goes to the soil mineral N pool.
! This could be updated later to divide the inputs between mineral N absorbed
! directly into the canopy and mineral N entering the soil pool.
!
! !USES:
   use clmtype
!  use clm_atmlnd   , only : clm_a2l
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
!
! !CALLED FROM:
! subroutine CNEcosystemDyn, in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 6/1/04: Created by Peter Thornton
! 11/06/09: Copy to all columns NOT just over soil. S. Levis
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   real, pointer :: forc_ndep(:)  ! nitrogen deposition rate (gN/m2/s)
   integer , pointer :: gridcell(:)   ! index into gridcell level quantities
!
! local pointers to implicit out scalars
!
   real, pointer :: ndep_to_sminn(:)
!
! !OTHER LOCAL VARIABLES:
   integer :: g,c                    ! indices

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
   forc_ndep     => clm3%g%forc_ndep    ! nitrogen deposition rate (gN/m2/s)
   gridcell      => clm3%g%l%c%gridcell

   ! Assign local pointers to derived type arrays (out)
   ndep_to_sminn => clm3%g%l%c%cnf%ndep_to_sminn

   ! Loop through columns
   do c = lbc, ubc
      g = gridcell(c)

      ndep_to_sminn(c) = forc_ndep(g)
      
   end do

end subroutine CNNDeposition

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNNFixation
!
! !INTERFACE:
subroutine CNNFixation(num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, update the nitrogen fixation rate
! as a function of annual total NPP. This rate gets updated once per year.
! All N fixation goes to the soil mineral N pool.
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine CNEcosystemDyn, in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 6/1/04: Created by Peter Thornton
! 2/14/05, PET: After looking at a number of point simulations,
!               it looks like a constant Nfix might be more efficient and 
!               maybe more realistic - setting to constant 0.4 gN/m2/yr.
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   real, pointer :: cannsum_npp(:) ! nitrogen deposition rate (gN/m2/s)
!
! local pointers to implicit out scalars
!
   real, pointer :: nfix_to_sminn(:)
!
! !OTHER LOCAL VARIABLES:
   integer  :: c,fc                  ! indices
   real :: t                     ! temporary

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
   cannsum_npp   => clm3%g%l%c%cps%cannsum_npp

   ! Assign local pointers to derived type arrays (out)
   nfix_to_sminn => clm3%g%l%c%cnf%nfix_to_sminn

   ! Loop through columns
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! the value 0.001666 is set to give 100 TgN/yr when global
      ! NPP = 60 PgC/yr.  (Cleveland et al., 1999)
      ! Convert from gN/m2/yr -> gN/m2/s
      !t = cannsum_npp(c) * 0.001666 / (86400. * 365.)
      t = (1.8 * (1. - exp(-0.003 * cannsum_npp(c))))/(86400. * 365.)
      nfix_to_sminn(c) = max(0.,t)
      ! PET 2/14/05: commenting out the dependence on NPP, and
      ! forcing Nfix to global constant = 0.4 gN/m2/yr
      !nfix_to_sminn(c) = 0.4 / (86400.*365.)

   end do

end subroutine CNNFixation

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNNLeaching
!
! !INTERFACE:
subroutine CNNLeaching(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, update the nitrogen leaching rate
! as a function of soluble mineral N and total soil water outflow.
!
! !USES:
   use clmtype
   use clm_varpar      , only : nlevsoi
   use clm_time_manager    , only : get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 6/9/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   real, pointer :: h2osoi_liq(:)    ! column liquid water (kg/m2) (new)
   real, pointer :: qflx_drain(:)    ! sub-surface runoff (mm H2O /s)
   real, pointer :: sminn(:)         ! (gN/m2) soil mineral N
!
! local pointers to implicit out scalars
!
   real, pointer :: sminn_leached(:) ! rate of mineral N leaching (gN/m2/s)
!
! !OTHER LOCAL VARIABLES:
   integer  :: j,c,fc             ! indices
   real :: dt                 ! radiation time step (seconds)
   real :: sf                 ! soluble fraction of mineral N (unitless)
   real :: disn_conc          ! dissolved mineral N concentration
                                  ! (gN/kg water)

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
   h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
   qflx_drain    => clm3%g%l%c%cwf%qflx_drain
   sminn         => clm3%g%l%c%cns%sminn

   ! Assign local pointers to derived type arrays (out)
   sminn_leached => clm3%g%l%c%cnf%sminn_leached

   ! set time steps
   dt = real( get_step_size() )

   ! Assume that 10% of the soil mineral N is in a soluble form
   sf = 0.1

   ! Loop through columns
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! calculate the dissolved mineral N concentration (gN/kg water)
      ! assumes that 10% of mineral nitrogen is soluble
      disn_conc = 0.
      if (h2osoi_liq(c) > 0.) then
         disn_conc = (sf * sminn(c))/h2osoi_liq(c)
      end if

      ! calculate the N leaching flux as a function of the dissolved
      ! concentration and the sub-surface drainage flux
      sminn_leached(c) = disn_conc * qflx_drain(c)

      ! limit the flux based on current sminn state
      ! only let at most the assumed soluble fraction
      ! of sminn be leached on any given timestep
      sminn_leached(c) = min(sminn_leached(c), (sf * sminn(c))/dt)
      
      ! limit the flux to a positive value
      sminn_leached(c) = max(sminn_leached(c), 0.)

   end do

end subroutine CNNLeaching

end module CNNDynamicsMod
