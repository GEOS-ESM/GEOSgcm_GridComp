module CNAnnualUpdateMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNAnnualUpdateMod
!
! !DESCRIPTION:
! Module for updating annual summation variables
!
! !USES:
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: CNAnnualUpdate
!
! !REVISION HISTORY:
! 4/23/2004: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNAnnualUpdate
!
! !INTERFACE:
subroutine CNAnnualUpdate(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                          num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update annual summation variables
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size, get_days_per_year
   use clm_varcon      , only: secspday
   use subgridAveMod   , only: p2c
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: lbp, ubp        ! pft bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(ubp-lbp+1) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 10/1/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: pcolumn(:)               ! index into column level
                                                 ! quantities
!
! local pointers to implicit in/out scalars
!
   real, pointer :: annsum_counter(:)        ! seconds since last annual accumulator turnover
   real, pointer :: tempsum_potential_gpp(:) ! temporary annual sum of potential GPP
   real, pointer :: annsum_potential_gpp(:)  ! annual sum of potential GPP
   real, pointer :: tempmax_retransn(:)      ! temporary annual max of retranslocated N pool (gN/m2)
   real, pointer :: annmax_retransn(:)       ! annual max of retranslocated N pool (gN/m2)
   real, pointer :: tempavg_t2m(:)           ! temporary average 2m air temperature (K)
   real, pointer :: annavg_t2m(:)            ! annual average 2m air temperature (K)
   real, pointer :: tempsum_npp(:)           ! temporary sum NPP (gC/m2/yr)
   real, pointer :: annsum_npp(:)            ! annual sum NPP (gC/m2/yr)
   real, pointer :: cannsum_npp(:)           ! column annual sum NPP (gC/m2/yr)
   real, pointer :: cannavg_t2m(:)    !annual average of 2m air temperature, averaged from pft-level (K)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p          ! indices
   integer :: fp,fc        ! lake filter indices
   real:: dt           ! radiation time step (seconds)

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays
   annsum_counter        => clm3%g%l%c%cps%annsum_counter
   tempsum_potential_gpp => clm3%g%l%c%p%pepv%tempsum_potential_gpp
   annsum_potential_gpp  => clm3%g%l%c%p%pepv%annsum_potential_gpp
   tempmax_retransn      => clm3%g%l%c%p%pepv%tempmax_retransn
   annmax_retransn       => clm3%g%l%c%p%pepv%annmax_retransn
   tempavg_t2m           => clm3%g%l%c%p%pepv%tempavg_t2m
   annavg_t2m            => clm3%g%l%c%p%pepv%annavg_t2m
   tempsum_npp           => clm3%g%l%c%p%pepv%tempsum_npp
   annsum_npp            => clm3%g%l%c%p%pepv%annsum_npp
   cannsum_npp           => clm3%g%l%c%cps%cannsum_npp
   cannavg_t2m           => clm3%g%l%c%cps%cannavg_t2m
   pcolumn               => clm3%g%l%c%p%column

   ! set time steps
   dt = real( get_step_size() )

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      annsum_counter(c) = annsum_counter(c) + dt
   end do

!!!if (annsum_counter(filter_soilc(1)) >= get_days_per_year() * secspday) then ! new (slevis)

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! In the future -- REMOVE this code and use the equivalent code above always
      c = pcolumn(p)                                                ! old (slevis)
      if (annsum_counter(c) >= get_days_per_year() * secspday) then ! old (slevis)

         ! update annual plant ndemand accumulator
         annsum_potential_gpp(p)  = tempsum_potential_gpp(p)
         tempsum_potential_gpp(p) = 0.

         ! update annual total N retranslocation accumulator
         annmax_retransn(p)  = tempmax_retransn(p)
         tempmax_retransn(p) = 0.

         ! update annual average 2m air temperature accumulator
         annavg_t2m(p)  = tempavg_t2m(p)
         tempavg_t2m(p) = 0.

         ! update annual NPP accumulator, convert to annual total
         annsum_npp(p) = tempsum_npp(p) * dt
         tempsum_npp(p) = 0.

      end if ! old (slevis)
   end do

   ! use p2c routine to get selected column-average pft-level fluxes and states
   call p2c(num_soilc, filter_soilc, annsum_npp, cannsum_npp)
   call p2c(num_soilc, filter_soilc, annavg_t2m, cannavg_t2m)

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      if (annsum_counter(c) >= get_days_per_year() * 86400.) annsum_counter(c) = 0.
   end do

end subroutine CNAnnualUpdate
!-----------------------------------------------------------------------

end module CNAnnualUpdateMod
