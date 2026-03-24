module CNPrecisionControlMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: CNPrecisionControlMod
! 
! !DESCRIPTION:
! controls on very low values in critical state variables 
! 
! !USES:
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: CNPrecisionControl
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
! !IROUTINE: CNPrecisionControl
!
! !INTERFACE:
subroutine CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION: 
! On the radiation time step, force leaf and deadstem c and n to 0 if
! they get too small.
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 8/1/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   real, pointer :: col_ctrunc(:)         ! (gC/m2) column-level sink for C truncation
   real, pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real, pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real, pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real, pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real, pointer :: soil1c(:)             ! (gC/m2) soil organic matter C (fast pool)
   real, pointer :: soil2c(:)             ! (gC/m2) soil organic matter C (medium pool)
   real, pointer :: soil3c(:)             ! (gC/m2) soil organic matter C (slow pool)
   real, pointer :: soil4c(:)             ! (gC/m2) soil organic matter C (slowest pool)
   real, pointer :: col_ntrunc(:)         ! (gN/m2) column-level sink for N truncation
   real, pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real, pointer :: litr1n(:)             ! (gN/m2) litter labile N
   real, pointer :: litr2n(:)             ! (gN/m2) litter cellulose N
   real, pointer :: litr3n(:)             ! (gN/m2) litter lignin N
   real, pointer :: soil1n(:)             ! (gN/m2) soil organic matter N (fast pool)
   real, pointer :: soil2n(:)             ! (gN/m2) soil organic matter N (medium pool)
   real, pointer :: soil3n(:)             ! (gN/m2) soil orgainc matter N (slow pool)
   real, pointer :: soil4n(:)             ! (gN/m2) soil orgainc matter N (slowest pool)
   real, pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real, pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real, pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real, pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real, pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real, pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real, pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real, pointer :: frootc(:)             ! (gC/m2) fine root C
   real, pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real, pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real, pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real, pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real, pointer :: leafc(:)              ! (gC/m2) leaf C
   real, pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real, pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real, pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real, pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real, pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real, pointer :: livestemc(:)          ! (gC/m2) live stem C
   real, pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real, pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real, pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
   real, pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real, pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real, pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real, pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real, pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real, pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real, pointer :: frootn(:)             ! (gN/m2) fine root N
   real, pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real, pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real, pointer :: leafn(:)              ! (gN/m2) leaf N 
   real, pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real, pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real, pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real, pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real, pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real, pointer :: livestemn(:)          ! (gN/m2) live stem N
   real, pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real, pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real, pointer :: npool(:)              ! (gN/m2) temporary plant N pool
   real, pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
   real, pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
!
! local pointers to implicit in/out scalars
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p      ! indices
   integer :: fp,fc    ! lake filter indices
   real:: pc,pn    ! truncation terms for pft-level corrections
   real:: cc,cn    ! truncation terms for column-level corrections
   real:: ccrit    ! critical carbon state value for truncation
   real:: ncrit    ! critical nitrogen state value for truncation
    
!EOP
!-----------------------------------------------------------------------
    ! assign local pointers at the column level
    col_ctrunc                     => clm3%g%l%c%ccs%col_ctrunc
    cwdc                           => clm3%g%l%c%ccs%cwdc
    litr1c                         => clm3%g%l%c%ccs%litr1c
    litr2c                         => clm3%g%l%c%ccs%litr2c
    litr3c                         => clm3%g%l%c%ccs%litr3c
    soil1c                         => clm3%g%l%c%ccs%soil1c
    soil2c                         => clm3%g%l%c%ccs%soil2c
    soil3c                         => clm3%g%l%c%ccs%soil3c
    soil4c                         => clm3%g%l%c%ccs%soil4c
    col_ntrunc                     => clm3%g%l%c%cns%col_ntrunc
    cwdn                           => clm3%g%l%c%cns%cwdn
    litr1n                         => clm3%g%l%c%cns%litr1n
    litr2n                         => clm3%g%l%c%cns%litr2n
    litr3n                         => clm3%g%l%c%cns%litr3n
    soil1n                         => clm3%g%l%c%cns%soil1n
    soil2n                         => clm3%g%l%c%cns%soil2n
    soil3n                         => clm3%g%l%c%cns%soil3n
    soil4n                         => clm3%g%l%c%cns%soil4n

    ! assign local pointers at the pft level
    cpool                          => clm3%g%l%c%p%pcs%cpool
    deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
    deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
    deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
    deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
    deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
    deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
    frootc                         => clm3%g%l%c%p%pcs%frootc
    frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
    frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
    gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
    gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
    leafc                          => clm3%g%l%c%p%pcs%leafc
    leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
    leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
    livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
    livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
    livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
    livestemc                      => clm3%g%l%c%p%pcs%livestemc
    livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
    livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
    pft_ctrunc                     => clm3%g%l%c%p%pcs%pft_ctrunc
    deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    deadstemn                      => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    frootn                         => clm3%g%l%c%p%pns%frootn
    frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    leafn                          => clm3%g%l%c%p%pns%leafn
    leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    livecrootn                     => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    livestemn                      => clm3%g%l%c%p%pns%livestemn
    livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    npool                          => clm3%g%l%c%p%pns%npool
    pft_ntrunc                     => clm3%g%l%c%p%pns%pft_ntrunc
    retransn                       => clm3%g%l%c%p%pns%retransn
   
   ! set the critical carbon state value for truncation (gC/m2)
   ccrit = 1.e-8
   ! set the critical nitrogen state value for truncation (gN/m2)
   ncrit = 1.e-8
   
   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      
      ! initialize the pft-level C and N truncation terms
      pc = 0.
      pn = 0.
      
      ! do tests on state variables for precision control
      ! for linked C-N state variables, perform precision test on
      ! the C component, but truncate C, C13, and N components
      
      ! leaf C and N
      if (abs(leafc(p)) < ccrit) then
          pc = pc + leafc(p)
          leafc(p) = 0.
          pn = pn + leafn(p)
          leafn(p) = 0.
      end if

      ! leaf storage C and N
      if (abs(leafc_storage(p)) < ccrit) then
          pc = pc + leafc_storage(p)
          leafc_storage(p) = 0.
          pn = pn + leafn_storage(p)
          leafn_storage(p) = 0.
      end if
          
      ! leaf transfer C and N
      if (abs(leafc_xfer(p)) < ccrit) then
          pc = pc + leafc_xfer(p)
          leafc_xfer(p) = 0.
          pn = pn + leafn_xfer(p)
          leafn_xfer(p) = 0.
      end if
          
      ! froot C and N
      if (abs(frootc(p)) < ccrit) then
          pc = pc + frootc(p)
          frootc(p) = 0.
          pn = pn + frootn(p)
          frootn(p) = 0.
      end if

      ! froot storage C and N
      if (abs(frootc_storage(p)) < ccrit) then
          pc = pc + frootc_storage(p)
          frootc_storage(p) = 0.
          pn = pn + frootn_storage(p)
          frootn_storage(p) = 0.
      end if
          
      ! froot transfer C and N
      if (abs(frootc_xfer(p)) < ccrit) then
          pc = pc + frootc_xfer(p)
          frootc_xfer(p) = 0.
          pn = pn + frootn_xfer(p)
          frootn_xfer(p) = 0.
      end if
          
      ! livestem C and N
      if (abs(livestemc(p)) < ccrit) then
          pc = pc + livestemc(p)
          livestemc(p) = 0.
          pn = pn + livestemn(p)
          livestemn(p) = 0.
      end if

      ! livestem storage C and N
      if (abs(livestemc_storage(p)) < ccrit) then
          pc = pc + livestemc_storage(p)
          livestemc_storage(p) = 0.
          pn = pn + livestemn_storage(p)
          livestemn_storage(p) = 0.
      end if
          
      ! livestem transfer C and N
      if (abs(livestemc_xfer(p)) < ccrit) then
          pc = pc + livestemc_xfer(p)
          livestemc_xfer(p) = 0.
          pn = pn + livestemn_xfer(p)
          livestemn_xfer(p) = 0.
      end if
          
      ! deadstem C and N
      if (abs(deadstemc(p)) < ccrit) then
          pc = pc + deadstemc(p)
          deadstemc(p) = 0.
          pn = pn + deadstemn(p)
          deadstemn(p) = 0.
      end if

      ! deadstem storage C and N
      if (abs(deadstemc_storage(p)) < ccrit) then
          pc = pc + deadstemc_storage(p)
          deadstemc_storage(p) = 0.
          pn = pn + deadstemn_storage(p)
          deadstemn_storage(p) = 0.
      end if
          
      ! deadstem transfer C and N
      if (abs(deadstemc_xfer(p)) < ccrit) then
          pc = pc + deadstemc_xfer(p)
          deadstemc_xfer(p) = 0.
          pn = pn + deadstemn_xfer(p)
          deadstemn_xfer(p) = 0.
      end if
          
      ! livecroot C and N
      if (abs(livecrootc(p)) < ccrit) then
          pc = pc + livecrootc(p)
          livecrootc(p) = 0.
          pn = pn + livecrootn(p)
          livecrootn(p) = 0.
      end if

      ! livecroot storage C and N
      if (abs(livecrootc_storage(p)) < ccrit) then
          pc = pc + livecrootc_storage(p)
          livecrootc_storage(p) = 0.
          pn = pn + livecrootn_storage(p)
          livecrootn_storage(p) = 0.
      end if
          
      ! livecroot transfer C and N
      if (abs(livecrootc_xfer(p)) < ccrit) then
          pc = pc + livecrootc_xfer(p)
          livecrootc_xfer(p) = 0.
          pn = pn + livecrootn_xfer(p)
          livecrootn_xfer(p) = 0.
      end if
          
      ! deadcroot C and N
      if (abs(deadcrootc(p)) < ccrit) then
          pc = pc + deadcrootc(p)
          deadcrootc(p) = 0.
          pn = pn + deadcrootn(p)
          deadcrootn(p) = 0.
      end if

      ! deadcroot storage C and N
      if (abs(deadcrootc_storage(p)) < ccrit) then
          pc = pc + deadcrootc_storage(p)
          deadcrootc_storage(p) = 0.
          pn = pn + deadcrootn_storage(p)
          deadcrootn_storage(p) = 0.
      end if
          
      ! deadcroot transfer C and N
      if (abs(deadcrootc_xfer(p)) < ccrit) then
          pc = pc + deadcrootc_xfer(p)
          deadcrootc_xfer(p) = 0.
          pn = pn + deadcrootn_xfer(p)
          deadcrootn_xfer(p) = 0.
      end if
          
      ! gresp_storage (C only)
      if (abs(gresp_storage(p)) < ccrit) then
          pc = pc + gresp_storage(p)
          gresp_storage(p) = 0.
      end if

      ! gresp_xfer (C only)
      if (abs(gresp_xfer(p)) < ccrit) then
          pc = pc + gresp_xfer(p)
          gresp_xfer(p) = 0.
      end if
          
      ! cpool (C only)
      if (abs(cpool(p)) < ccrit) then
          pc = pc + cpool(p)
          cpool(p) = 0.
      end if
          
      ! retransn (N only)
      if (abs(retransn(p)) < ncrit) then
          pn = pn + retransn(p)
          retransn(p) = 0.
      end if
          
      ! npool (N only)
      if (abs(npool(p)) < ncrit) then
          pn = pn + npool(p)
          npool(p) = 0.
      end if
      
      pft_ctrunc(p) = pft_ctrunc(p) + pc
      pft_ntrunc(p) = pft_ntrunc(p) + pn
          
   end do ! end of pft loop

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      
      ! initialize the column-level C and N truncation terms
      cc = 0.
      cn = 0.
      
      ! do tests on state variables for precision control
      ! for linked C-N state variables, perform precision test on
      ! the C component, but truncate both C and N components
      
      ! coarse woody debris C and N
      if (abs(cwdc(c)) < ccrit) then
          cc = cc + cwdc(c)
          cwdc(c) = 0.
          cn = cn + cwdn(c)
          cwdn(c) = 0.
      end if

      ! litr1 C and N
      if (abs(litr1c(c)) < ccrit) then
          cc = cc + litr1c(c)
          litr1c(c) = 0.
          cn = cn + litr1n(c)
          litr1n(c) = 0.
      end if

      ! litr2 C and N
      if (abs(litr2c(c)) < ccrit) then
          cc = cc + litr2c(c)
          litr2c(c) = 0.
          cn = cn + litr2n(c)
          litr2n(c) = 0.
      end if

      ! litr3 C and N
      if (abs(litr3c(c)) < ccrit) then
          cc = cc + litr3c(c)
          litr3c(c) = 0.
          cn = cn + litr3n(c)
          litr3n(c) = 0.
      end if

      ! soil1 C and N
      if (abs(soil1c(c)) < ccrit) then
          cc = cc + soil1c(c)
          soil1c(c) = 0.
          cn = cn + soil1n(c)
          soil1n(c) = 0.
      end if

      ! soil2 C and N
      if (abs(soil2c(c)) < ccrit) then
          cc = cc + soil2c(c)
          soil2c(c) = 0.
          cn = cn + soil2n(c)
          soil2n(c) = 0.
      end if

      ! soil3 C and N
      if (abs(soil3c(c)) < ccrit) then
          cc = cc + soil3c(c)
          soil3c(c) = 0.
          cn = cn + soil3n(c)
          soil3n(c) = 0.
      end if
      
      ! soil4 C and N
      if (abs(soil4c(c)) < ccrit) then
          cc = cc + soil4c(c)
          soil4c(c) = 0.
          cn = cn + soil4n(c)
          soil4n(c) = 0.
      end if
      
      ! not doing precision control on soil mineral N, since it will
      ! be getting the N truncation flux anyway.
      
      col_ctrunc(c) = col_ctrunc(c) + cc
      col_ntrunc(c) = col_ntrunc(c) + cn
      
   end do   ! end of column loop

end subroutine CNPrecisionControl
!-----------------------------------------------------------------------

end module CNPrecisionControlMod
