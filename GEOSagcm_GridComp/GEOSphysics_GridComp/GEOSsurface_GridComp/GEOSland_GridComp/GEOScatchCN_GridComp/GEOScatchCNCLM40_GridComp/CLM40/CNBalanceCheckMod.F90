module CNBalanceCheckMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNBalanceCheckMod
!
! !DESCRIPTION:
! Module for carbon mass balance checking.
!
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: BeginCBalance
    public :: BeginNBalance
    public :: CBalanceCheck
    public :: NBalanceCheck
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
! !IROUTINE: BeginCBalance
!
! !INTERFACE:
subroutine BeginCBalance(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, calculate the beginning carbon balance for mass
! conservation checks.
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 2/4/05: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real, pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
!
! local pointers to implicit out arrays
   real, pointer :: col_begcb(:)   ! carbon mass, beginning of time step (gC/m**2)
!
! !OTHER LOCAL VARIABLES:
   integer :: c     ! indices
   integer :: fc   ! lake filter indices
!
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers at the column level
   col_begcb                      => clm3%g%l%c%ccbal%begcb
   totcolc                        => clm3%g%l%c%ccs%totcolc

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
 
      ! calculate beginning column-level carbon balance,
      ! for mass conservation check
 
      col_begcb(c) = totcolc(c)

   end do ! end of columns loop
 

end subroutine BeginCBalance
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BeginNBalance
!
! !INTERFACE:
subroutine BeginNBalance(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, calculate the beginning nitrogen balance for mass
! conservation checks.
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 2/4/05: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real, pointer :: totcoln(:)            ! (gN/m2) total column nitrogen, incl veg
!
! local pointers to implicit out arrays
   real, pointer :: col_begnb(:)   ! nitrogen mass, beginning of time step (gN/m**2)
!
! !OTHER LOCAL VARIABLES:
   integer :: c     ! indices
   integer :: fc   ! lake filter indices
!
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers at the column level
   col_begnb                      => clm3%g%l%c%cnbal%begnb
   totcoln                        => clm3%g%l%c%cns%totcoln

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
 
      ! calculate beginning column-level nitrogen balance,
      ! for mass conservation check
 
      col_begnb(c) = totcoln(c)

   end do ! end of columns loop
 
end subroutine BeginNBalance
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CBalanceCheck
!
! !INTERFACE:
subroutine CBalanceCheck(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, perform carbon mass conservation check for column and pft
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 12/9/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arrays
   real, pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
   real, pointer :: gpp(:)            ! (gC/m2/s) gross primary production 
   real, pointer :: er(:)            ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
   real, pointer :: col_fire_closs(:) ! (gC/m2/s) total column-level fire C loss
   real, pointer :: col_hrv_xsmrpool_to_atm(:)  ! excess MR pool harvest mortality (gC/m2/s)
   real, pointer :: dwt_closs(:)              ! (gC/m2/s) total carbon loss from product pools and conversion
   real, pointer :: product_closs(:)      ! (gC/m2/s) total wood product carbon loss
!
! local pointers to implicit out arrays
   real, pointer :: col_cinputs(:)  ! (gC/m2/s) total column-level carbon inputs (for balance check)
   real, pointer :: col_coutputs(:) ! (gC/m2/s) total column-level carbon outputs (for balance check)
   real, pointer :: col_begcb(:)    ! carbon mass, beginning of time step (gC/m**2)
   real, pointer :: col_endcb(:)    ! carbon mass, end of time step (gC/m**2)
   real, pointer :: col_errcb(:)    ! carbon balance error for the timestep (gC/m**2)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,err_index  ! indices
   integer :: fc          ! lake filter indices
   logical :: err_found      ! error flag
   real:: dt             ! radiation time step (seconds)
   integer :: icnt           ! counter
!EOP
!-----------------------------------------------------------------------

    ! assign local pointers to column-level arrays
	totcolc                        => clm3%g%l%c%ccs%totcolc
	gpp                            => clm3%g%l%c%ccf%pcf_a%gpp
	er                             => clm3%g%l%c%ccf%er
	col_fire_closs                 => clm3%g%l%c%ccf%col_fire_closs
	col_hrv_xsmrpool_to_atm        => clm3%g%l%c%ccf%pcf_a%hrv_xsmrpool_to_atm
	dwt_closs                      => clm3%g%l%c%ccf%dwt_closs
	product_closs                  => clm3%g%l%c%ccf%product_closs
	
    col_cinputs                    => clm3%g%l%c%ccf%col_cinputs
    col_coutputs                   => clm3%g%l%c%ccf%col_coutputs
    col_begcb                      => clm3%g%l%c%ccbal%begcb
    col_endcb                      => clm3%g%l%c%ccbal%endcb
    col_errcb                      => clm3%g%l%c%ccbal%errcb

   ! set time steps
   dt = real( get_step_size() )

   icnt = 0
   err_found = .false.
   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! calculate the total column-level carbon storage, for mass conservation check

      col_endcb(c) = totcolc(c)

      ! calculate total column-level inputs
	  
	  col_cinputs(c) = gpp(c)
	  
      ! calculate total column-level outputs
	  ! er = ar + hr, col_fire_closs includes pft-level fire losses

      col_coutputs(c) = er(c) + col_fire_closs(c) + dwt_closs(c) + product_closs(c) + col_hrv_xsmrpool_to_atm(c)

      ! calculate the total column-level carbon balance error for this time step

      col_errcb(c) = (col_cinputs(c) - col_coutputs(c))*dt - &
         (col_endcb(c) - col_begcb(c))

      ! check for significant errors (1e-8 for real*8, 1e0 for real*4)
      if (col_endcb(c) > 0. .and. abs(col_errcb(c)) > 0.1) then
         err_found = .true.
         err_index = c
         icnt = icnt + 1
         if(icnt > 1 .and. abs(col_errcb(c)) > abs(col_errcb(err_index))) err_index = c
      endif
      
   end do ! end of columns loop

   if (err_found) then
      c = err_index
      write(99,*)'column cbalance error = ', col_errcb(c), c, icnt
      write(99,*)'begcb       = ',col_begcb(c)
      write(99,*)'endcb       = ',col_endcb(c)
      write(99,*)'delta store = ',col_endcb(c)-col_begcb(c)
      write(99,*)'input mass  = ',col_cinputs(c)*dt
      write(99,*)'output mass = ',col_coutputs(c)*dt
      write(99,*)'net flux    = ',(col_cinputs(c)-col_coutputs(c))*dt
	  write(99,*)'nee         = ',clm3%g%l%c%ccf%nee(c) * dt
	  write(99,*)'gpp         = ',gpp(c) * dt
	  write(99,*)'er          = ',er(c) * dt
	  write(99,*)'col_fire_closs         = ',col_fire_closs(c) * dt
     write(99,*)'col_hrv_xsmrpool_to_atm = ',col_hrv_xsmrpool_to_atm(c) * dt
	  write(99,*)'dwt_closs         = ',dwt_closs(c) * dt
	  write(99,*)'product_closs         = ',product_closs(c) * dt
!!    stop 'C balance' ! call endrun ! gkw
   end if


end subroutine CBalanceCheck
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NBalanceCheck
!
! !INTERFACE:
subroutine NBalanceCheck(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, perform nitrogen mass conservation check
! for column and pft
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 12/9/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arrays
   real, pointer :: totcoln(:)               ! (gN/m2) total column nitrogen, incl veg
   real, pointer :: ndep_to_sminn(:)         ! atmospheric N deposition to soil mineral N (gN/m2/s)
   real, pointer :: nfix_to_sminn(:)         ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s) 
   real, pointer :: supplement_to_sminn(:)   ! supplemental N supply (gN/m2/s)
   real, pointer :: denit(:)                 ! total rate of denitrification (gN/m2/s)
   real, pointer :: sminn_leached(:)         ! soil mineral N pool loss to leaching (gN/m2/s)
   real, pointer :: col_fire_nloss(:)        ! total column-level fire N loss (gN/m2/s)
   real, pointer :: dwt_nloss(:)              ! (gN/m2/s) total nitrogen loss from product pools and conversion
   real, pointer :: product_nloss(:)          ! (gN/m2/s) total wood product nitrogen loss
!
! local pointers to implicit in/out arrays
!
! local pointers to implicit out arrays
   real, pointer :: col_ninputs(:)           ! column-level N inputs (gN/m2/s)
   real, pointer :: col_noutputs(:)          ! column-level N outputs (gN/m2/s)
   real, pointer :: col_begnb(:)             ! nitrogen mass, beginning of time step (gN/m**2)
   real, pointer :: col_endnb(:)             ! nitrogen mass, end of time step (gN/m**2)
   real, pointer :: col_errnb(:)             ! nitrogen balance error for the timestep (gN/m**2)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,err_index    ! indices
   integer :: fc             ! lake filter indices
   logical :: err_found      ! error flag
   real:: dt             ! radiation time step (seconds)
   integer :: icnt           ! counter
!EOP
!-----------------------------------------------------------------------
    ! assign local pointers to column-level arrays

    totcoln                        => clm3%g%l%c%cns%totcoln
	ndep_to_sminn                  => clm3%g%l%c%cnf%ndep_to_sminn
    nfix_to_sminn                  => clm3%g%l%c%cnf%nfix_to_sminn
    supplement_to_sminn            => clm3%g%l%c%cnf%supplement_to_sminn
    denit                          => clm3%g%l%c%cnf%denit
    sminn_leached                  => clm3%g%l%c%cnf%sminn_leached
    col_fire_nloss                 => clm3%g%l%c%cnf%col_fire_nloss
    dwt_nloss                      => clm3%g%l%c%cnf%dwt_nloss
    product_nloss                  => clm3%g%l%c%cnf%product_nloss

    col_ninputs                    => clm3%g%l%c%cnf%col_ninputs
    col_noutputs                   => clm3%g%l%c%cnf%col_noutputs
    col_begnb                      => clm3%g%l%c%cnbal%begnb
    col_endnb                      => clm3%g%l%c%cnbal%endnb
    col_errnb                      => clm3%g%l%c%cnbal%errnb

   ! set time steps
   dt = real( get_step_size() )

   err_found = .false.
   ! column loop
   icnt = 0
   do fc = 1,num_soilc
      c=filter_soilc(fc)

      ! calculate the total column-level nitrogen storage, for mass conservation check

      col_endnb(c) = totcoln(c)

      ! calculate total column-level inputs

      col_ninputs(c) = ndep_to_sminn(c) + nfix_to_sminn(c) + supplement_to_sminn(c)

      ! calculate total column-level outputs

      col_noutputs(c) = denit(c) + sminn_leached(c) + col_fire_nloss(c) + dwt_nloss(c) + product_nloss(c)

      ! calculate the total column-level nitrogen balance error for this time step

      col_errnb(c) = (col_ninputs(c) - col_noutputs(c))*dt - &
         (col_endnb(c) - col_begnb(c))

      if (col_endnb(c) > 0. .and. abs(col_errnb(c)) > 0.1) then
         err_found = .true.
         err_index = c
         icnt = icnt + 1
         if(icnt > 1 .and. abs(col_errnb(c)) > abs(col_errnb(err_index))) err_index = c
      endif

   end do ! end of columns loop

   if (err_found) then
      c = err_index
      write(99,*)'column nbalance error = ', col_errnb(c), c
      write(99,*)'begnb       = ',col_begnb(c)
      write(99,*)'endnb       = ',col_endnb(c)
      write(99,*)'delta store = ',col_endnb(c)-col_begnb(c)
      write(99,*)'input mass  = ',col_ninputs(c)*dt
      write(99,*)'output mass = ',col_noutputs(c)*dt
      write(99,*)'net flux    = ',(col_ninputs(c)-col_noutputs(c))*dt
!!    stop 'N balance' ! call endrun ! gkw
   end if

end subroutine NBalanceCheck
!-----------------------------------------------------------------------

end module CNBalanceCheckMod
