module CNAllocationMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNAllocationMod
!
! !DESCRIPTION:
! Module holding routines used in allocation model for coupled carbon
! nitrogen code.
!
! !USES:
  implicit none
  save
  private
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNAllocation
!
! !REVISION HISTORY:
! 8/5/03: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNAllocation
!
! !INTERFACE:
subroutine CNAllocation (lbp, ubp, lbc, ubc, &
       num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use subgridAveMod, only: p2c
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbp, ubp        ! pft-index bounds
   integer, intent(in) :: lbc, ubc        ! column-index bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNdecompAlloc in module CNdecompMod.F90
!
! !REVISION HISTORY:
! 8/5/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: migrated to vector data structures
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   ! pft level
   integer , pointer :: ivt(:)        ! pft vegetation type
   integer , pointer :: pcolumn(:)    ! pft's column index
   real, pointer :: lgsf(:)       ! long growing season factor [0-1]
   real, pointer :: xsmrpool(:)      ! (kgC/m2) temporary photosynthate C pool
   real, pointer :: retransn(:)   ! (kgN/m2) plant pool of retranslocated N
   real, pointer :: psnsun(:)     ! sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
   real, pointer :: psnsha(:)     ! shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
   real, pointer :: laisun(:)     ! sunlit projected leaf area index
   real, pointer :: laisha(:)     ! shaded projected leaf area index
   real, pointer :: leaf_mr(:)
   real, pointer :: froot_mr(:)
   real, pointer :: livestem_mr(:)
   real, pointer :: livecroot_mr(:)
   real, pointer :: leaf_curmr(:)
   real, pointer :: froot_curmr(:)
   real, pointer :: livestem_curmr(:)
   real, pointer :: livecroot_curmr(:)
   real, pointer :: leaf_xsmr(:)
   real, pointer :: froot_xsmr(:)
   real, pointer :: livestem_xsmr(:)
   real, pointer :: livecroot_xsmr(:)
   ! column level
   real, pointer :: sminn(:)      ! (kgN/m2) soil mineral N
   ! ecophysiological constants
   real, pointer :: woody(:)      ! binary flag for woody lifeform (1=woody, 0=not woody)
   real, pointer :: froot_leaf(:) ! allocation parameter: new fine root C per new leaf C (gC/gC)
   real, pointer :: croot_stem(:) ! allocation parameter: new coarse root C per new stem C (gC/gC)
   real, pointer :: stem_leaf(:)  ! allocation parameter: new stem c per new leaf C (gC/gC)
   real, pointer :: flivewd(:)    ! allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
   real, pointer :: leafcn(:)     ! leaf C:N (gC/gN)
   real, pointer :: frootcn(:)    ! fine root C:N (gC/gN)
   real, pointer :: livewdcn(:)   ! live wood (phloem and ray parenchyma) C:N (gC/gN)
   real, pointer :: deadwdcn(:)   ! dead wood (xylem and heartwood) C:N (gC/gN)
   real, pointer :: fcur2(:)      ! allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
   integer, pointer :: plandunit(:)   ! index into landunit level quantities
   integer, pointer :: clandunit(:)   ! index into landunit level quantities
   integer , pointer :: itypelun(:)   ! landunit type
!
! local pointers to implicit in/out arrays
!
   ! pft level
   real, pointer :: gpp(:)                   ! GPP flux before downregulation (gC/m2/s)
   real, pointer :: availc(:)                ! C flux available for allocation (gC/m2/s)
   real, pointer :: xsmrpool_recover(:)         ! C flux assigned to recovery of negative cpool (gC/m2/s)
   real, pointer :: c_allometry(:)           ! C allocation index (DIM)
   real, pointer :: n_allometry(:)           ! N allocation index (DIM)
   real, pointer :: plant_ndemand(:)         ! N flux required to support initial GPP (gN/m2/s)
   real, pointer :: tempsum_potential_gpp(:) ! temporary annual sum of potential GPP 
   real, pointer :: tempmax_retransn(:)      ! temporary annual max of retranslocated N pool (gN/m2)
   real, pointer :: annsum_potential_gpp(:)  ! annual sum of potential GPP
   real, pointer :: avail_retransn(:)        ! N flux available from retranslocation pool (gN/m2/s)
   real, pointer :: annmax_retransn(:)       ! annual max of retranslocated N pool
   real, pointer :: plant_nalloc(:)          ! total allocated N flux (gN/m2/s)
   real, pointer :: plant_calloc(:)          ! total allocated C flux (gC/m2/s)
   real, pointer :: excess_cflux(:)          ! C flux not allocated due to downregulation (gC/m2/s)
   real, pointer :: downreg(:)               ! fractional reduction in GPP due to N limitation (DIM)
   real, pointer :: annsum_npp(:)            ! annual sum of NPP, for wood allocation
   real, pointer :: cpool_to_xsmrpool(:)
   real, pointer :: psnsun_to_cpool(:)
   real, pointer :: psnshade_to_cpool(:)
   real, pointer :: cpool_to_leafc(:)
   real, pointer :: cpool_to_leafc_storage(:)
   real, pointer :: cpool_to_frootc(:)
   real, pointer :: cpool_to_frootc_storage(:)
   real, pointer :: cpool_to_livestemc(:)
   real, pointer :: cpool_to_livestemc_storage(:)
   real, pointer :: cpool_to_deadstemc(:)
   real, pointer :: cpool_to_deadstemc_storage(:)
   real, pointer :: cpool_to_livecrootc(:)
   real, pointer :: cpool_to_livecrootc_storage(:)
   real, pointer :: cpool_to_deadcrootc(:)
   real, pointer :: cpool_to_deadcrootc_storage(:)
   real, pointer :: cpool_to_gresp_storage(:)
   real, pointer :: retransn_to_npool(:)
   real, pointer :: sminn_to_npool(:)
   real, pointer :: npool_to_leafn(:)
   real, pointer :: npool_to_leafn_storage(:)
   real, pointer :: npool_to_frootn(:)
   real, pointer :: npool_to_frootn_storage(:)
   real, pointer :: npool_to_livestemn(:)
   real, pointer :: npool_to_livestemn_storage(:)
   real, pointer :: npool_to_deadstemn(:)
   real, pointer :: npool_to_deadstemn_storage(:)
   real, pointer :: npool_to_livecrootn(:)
   real, pointer :: npool_to_livecrootn_storage(:)
   real, pointer :: npool_to_deadcrootn(:)
   real, pointer :: npool_to_deadcrootn_storage(:)
   ! column level
   real, pointer :: fpi(:) ! fraction of potential immobilization (no units)
   real, pointer :: fpg(:) ! fraction of potential gpp (no units)
   real, pointer :: potential_immob(:)
   real, pointer :: actual_immob(:)
   real, pointer :: sminn_to_plant(:)
   real, pointer :: sminn_to_denit_excess(:)
   real, pointer :: supplement_to_sminn(:)
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p                  !indices
   integer :: fp                   !lake filter pft index
   integer :: fc                   !lake filter column index
   real:: dt                   !decomp timestep (seconds)
   integer :: nlimit               !flag for N limitation
   real, pointer:: col_plant_ndemand(:)    !column-level plant N demand
   real:: dayscrecover         !number of days to recover negative cpool
   real:: mr                   !maintenance respiration (gC/m2/s)
   real:: f1,f2,f3,f4,g1,g2    !allocation parameters
   real:: cnl,cnfr,cnlw,cndw   !C:N ratios for leaf, fine root, and wood
   real:: grperc, grpnow       !growth respirarion parameters
   real:: fcur                 !fraction of current psn displayed as growth
   real:: sum_ndemand          !total column N demand (gN/m2/s)
   real:: gresp_storage        !temporary variable for growth resp to storage
   real:: nlc                  !temporary variable for total new leaf carbon allocation
   real:: bdnr                 !bulk denitrification rate (1/s)
   real:: curmr, curmr_ratio   !xsmrpool temporary variables

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
   ivt                         => clm3%g%l%c%p%itype
   pcolumn                     => clm3%g%l%c%p%column
   plandunit                   => clm3%g%l%c%p%landunit
   clandunit                   => clm3%g%l%c%landunit
   itypelun                    => clm3%g%l%itype
   lgsf                        => clm3%g%l%c%p%pepv%lgsf
   xsmrpool                    => clm3%g%l%c%p%pcs%xsmrpool
   retransn                    => clm3%g%l%c%p%pns%retransn
   psnsun                      => clm3%g%l%c%p%pcf%psnsun
   psnsha                      => clm3%g%l%c%p%pcf%psnsha
   laisun                      => clm3%g%l%c%p%pps%laisun
   laisha                      => clm3%g%l%c%p%pps%laisha
   leaf_mr                     => clm3%g%l%c%p%pcf%leaf_mr
   froot_mr                    => clm3%g%l%c%p%pcf%froot_mr
   livestem_mr                 => clm3%g%l%c%p%pcf%livestem_mr
   livecroot_mr                => clm3%g%l%c%p%pcf%livecroot_mr
   leaf_curmr                  => clm3%g%l%c%p%pcf%leaf_curmr
   froot_curmr                 => clm3%g%l%c%p%pcf%froot_curmr
   livestem_curmr              => clm3%g%l%c%p%pcf%livestem_curmr
   livecroot_curmr             => clm3%g%l%c%p%pcf%livecroot_curmr
   leaf_xsmr                   => clm3%g%l%c%p%pcf%leaf_xsmr
   froot_xsmr                  => clm3%g%l%c%p%pcf%froot_xsmr
   livestem_xsmr               => clm3%g%l%c%p%pcf%livestem_xsmr
   livecroot_xsmr              => clm3%g%l%c%p%pcf%livecroot_xsmr
   sminn                       => clm3%g%l%c%cns%sminn
   woody                       => pftcon%woody
   froot_leaf                  => pftcon%froot_leaf
   croot_stem                  => pftcon%croot_stem
   stem_leaf                   => pftcon%stem_leaf
   flivewd                     => pftcon%flivewd
   leafcn                      => pftcon%leafcn
   frootcn                     => pftcon%frootcn
   livewdcn                    => pftcon%livewdcn
   deadwdcn                    => pftcon%deadwdcn
   fcur2                       => pftcon%fcur
   ! Assign local pointers to derived type arrays (out)
   gpp                         => clm3%g%l%c%p%pepv%gpp
   availc                      => clm3%g%l%c%p%pepv%availc
   xsmrpool_recover            => clm3%g%l%c%p%pepv%xsmrpool_recover
   c_allometry                 => clm3%g%l%c%p%pepv%c_allometry
   n_allometry                 => clm3%g%l%c%p%pepv%n_allometry
   plant_ndemand               => clm3%g%l%c%p%pepv%plant_ndemand
   tempsum_potential_gpp       => clm3%g%l%c%p%pepv%tempsum_potential_gpp
   tempmax_retransn            => clm3%g%l%c%p%pepv%tempmax_retransn
   annsum_potential_gpp        => clm3%g%l%c%p%pepv%annsum_potential_gpp
   avail_retransn              => clm3%g%l%c%p%pepv%avail_retransn
   annmax_retransn             => clm3%g%l%c%p%pepv%annmax_retransn
   plant_nalloc                => clm3%g%l%c%p%pepv%plant_nalloc
   plant_calloc                => clm3%g%l%c%p%pepv%plant_calloc
   excess_cflux                => clm3%g%l%c%p%pepv%excess_cflux
   downreg                     => clm3%g%l%c%p%pepv%downreg
   annsum_npp                  => clm3%g%l%c%p%pepv%annsum_npp
   cpool_to_xsmrpool           => clm3%g%l%c%p%pcf%cpool_to_xsmrpool
   psnsun_to_cpool             => clm3%g%l%c%p%pcf%psnsun_to_cpool
   psnshade_to_cpool           => clm3%g%l%c%p%pcf%psnshade_to_cpool
   cpool_to_leafc              => clm3%g%l%c%p%pcf%cpool_to_leafc
   cpool_to_leafc_storage      => clm3%g%l%c%p%pcf%cpool_to_leafc_storage
   cpool_to_frootc             => clm3%g%l%c%p%pcf%cpool_to_frootc
   cpool_to_frootc_storage     => clm3%g%l%c%p%pcf%cpool_to_frootc_storage
   cpool_to_livestemc          => clm3%g%l%c%p%pcf%cpool_to_livestemc
   cpool_to_livestemc_storage  => clm3%g%l%c%p%pcf%cpool_to_livestemc_storage
   cpool_to_deadstemc          => clm3%g%l%c%p%pcf%cpool_to_deadstemc
   cpool_to_deadstemc_storage  => clm3%g%l%c%p%pcf%cpool_to_deadstemc_storage
   cpool_to_livecrootc         => clm3%g%l%c%p%pcf%cpool_to_livecrootc
   cpool_to_livecrootc_storage => clm3%g%l%c%p%pcf%cpool_to_livecrootc_storage
   cpool_to_deadcrootc         => clm3%g%l%c%p%pcf%cpool_to_deadcrootc
   cpool_to_deadcrootc_storage => clm3%g%l%c%p%pcf%cpool_to_deadcrootc_storage
   cpool_to_gresp_storage      => clm3%g%l%c%p%pcf%cpool_to_gresp_storage
   retransn_to_npool           => clm3%g%l%c%p%pnf%retransn_to_npool
   sminn_to_npool              => clm3%g%l%c%p%pnf%sminn_to_npool
   npool_to_leafn              => clm3%g%l%c%p%pnf%npool_to_leafn
   npool_to_leafn_storage      => clm3%g%l%c%p%pnf%npool_to_leafn_storage
   npool_to_frootn             => clm3%g%l%c%p%pnf%npool_to_frootn
   npool_to_frootn_storage     => clm3%g%l%c%p%pnf%npool_to_frootn_storage
   npool_to_livestemn          => clm3%g%l%c%p%pnf%npool_to_livestemn
   npool_to_livestemn_storage  => clm3%g%l%c%p%pnf%npool_to_livestemn_storage
   npool_to_deadstemn          => clm3%g%l%c%p%pnf%npool_to_deadstemn
   npool_to_deadstemn_storage  => clm3%g%l%c%p%pnf%npool_to_deadstemn_storage
   npool_to_livecrootn         => clm3%g%l%c%p%pnf%npool_to_livecrootn
   npool_to_livecrootn_storage => clm3%g%l%c%p%pnf%npool_to_livecrootn_storage
   npool_to_deadcrootn         => clm3%g%l%c%p%pnf%npool_to_deadcrootn
   npool_to_deadcrootn_storage => clm3%g%l%c%p%pnf%npool_to_deadcrootn_storage
   fpi                         => clm3%g%l%c%cps%fpi
   fpg                         => clm3%g%l%c%cps%fpg
   potential_immob             => clm3%g%l%c%cnf%potential_immob
   actual_immob                => clm3%g%l%c%cnf%actual_immob
   sminn_to_plant              => clm3%g%l%c%cnf%sminn_to_plant
   sminn_to_denit_excess       => clm3%g%l%c%cnf%sminn_to_denit_excess
   supplement_to_sminn         => clm3%g%l%c%cnf%supplement_to_sminn

   ! set time steps
   dt = real( get_step_size() )

   ! set some space-and-time constant parameters 
   dayscrecover = 30.0
   grperc = 0.3
   grpnow = 1.0
   bdnr = 0.5 * (dt/86400.)

   ! loop over pfts to assess the total plant N demand
   do fp=1,num_soilp
      p = filter_soilp(fp)

      ! get the time step total gross photosynthesis
      ! this is coming from the canopy fluxes code, and is the
      ! gpp that is used to control stomatal conductance.
      ! For the nitrogen downregulation code, this is assumed
      ! to be the potential gpp, and the actual gpp will be
      ! reduced due to N limitation. 
      
      ! Convert psn from umol/m2/s -> gC/m2/s

      ! The input psn (psnsun and psnsha) are expressed per unit LAI
      ! in the sunlit and shaded canopy, respectively. These need to be
      ! scaled by laisun and laisha to get the total gpp for allocation

      psnsun_to_cpool(p) = psnsun(p) * laisun(p) * 12.011e-6
      psnshade_to_cpool(p) = psnsha(p) * laisha(p) * 12.011e-6
      
      gpp(p) = psnsun_to_cpool(p) + psnshade_to_cpool(p)

      ! get the time step total maintenance respiration
      ! These fluxes should already be in gC/m2/s

      mr = leaf_mr(p) + froot_mr(p)
      if (woody(ivt(p)) == 1.0) then
         mr = mr + livestem_mr(p) + livecroot_mr(p)
      end if

      ! carbon flux available for allocation
      availc(p) = gpp(p) - mr
      
      ! new code added for isotope calculations, 7/1/05, PET
      ! If mr > gpp, then some mr comes from gpp, the rest comes from
      ! cpool (xsmr)
      curmr_ratio = 1.
      if (mr > 0. .and. availc(p) < 0.) then
         curmr = gpp(p)
         curmr_ratio = curmr / mr
      end if
      leaf_curmr(p) = leaf_mr(p) * curmr_ratio
      leaf_xsmr(p) = leaf_mr(p) - leaf_curmr(p)
      froot_curmr(p) = froot_mr(p) * curmr_ratio
      froot_xsmr(p) = froot_mr(p) - froot_curmr(p)
      livestem_curmr(p) = livestem_mr(p) * curmr_ratio
      livestem_xsmr(p) = livestem_mr(p) - livestem_curmr(p)
      livecroot_curmr(p) = livecroot_mr(p) * curmr_ratio
      livecroot_xsmr(p) = livecroot_mr(p) - livecroot_curmr(p)
      
      ! no allocation when available c is negative
      availc(p) = max(availc(p),0.0)

      ! test for an xsmrpool deficit
      if (xsmrpool(p) < 0.0) then
         ! Running a deficit in the xsmrpool, so the first priority is to let
         ! some availc from this timestep accumulate in xsmrpool.
         ! Determine rate of recovery for xsmrpool deficit

         xsmrpool_recover(p) = -xsmrpool(p)/(dayscrecover*86400.0)
         if (xsmrpool_recover(p) < availc(p)) then
             ! available carbon reduced by amount for xsmrpool recovery
             availc(p) = availc(p) - xsmrpool_recover(p)
         else
             ! all of the available carbon goes to xsmrpool recovery
             xsmrpool_recover(p) = availc(p)
             availc(p) = 0.0
         end if
         cpool_to_xsmrpool(p) = xsmrpool_recover(p)
      end if

      f1 = froot_leaf(ivt(p))
      f2 = croot_stem(ivt(p))
     
      ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
      ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
      ! This variable allocation is only for trees. Shrubs have a constant
      ! allocation as specified in the pft-physiology file.  The value is also used
      ! as a trigger here: -1.0 means to use the dynamic allocation (trees).
      if (stem_leaf(ivt(p)) == -1.) then
         f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
      else
         f3 = stem_leaf(ivt(p))
      end if
      
      f4 = flivewd(ivt(p))
      g1 = grperc
      g2 = grpnow
      cnl = leafcn(ivt(p))
      cnfr = frootcn(ivt(p))
      cnlw = livewdcn(ivt(p))
      cndw = deadwdcn(ivt(p))

      ! based on available C, use constant allometric relationships to
      ! determine N requirements
      if (woody(ivt(p)) == 1.0) then
         c_allometry(p) = (1.+g1)*(1.+f1+f3*(1.+f2))
         n_allometry(p) = 1./cnl + f1/cnfr + (f3*f4*(1.+f2))/cnlw + &
                       (f3*(1.-f4)*(1.+f2))/cndw
      else
         c_allometry(p) = 1.+g1+f1+f1*g1
         n_allometry(p) = 1./cnl + f1/cnfr
      end if
      plant_ndemand(p) = availc(p)*(n_allometry(p)/c_allometry(p))

      ! retranslocated N deployment depends on seasonal cycle of potential GPP
      ! (requires one year run to accumulate demand)

      tempsum_potential_gpp(p) = tempsum_potential_gpp(p) + gpp(p)

      ! Adding the following line to carry max retransn info to CN Annual Update
      tempmax_retransn(p) = max(tempmax_retransn(p),retransn(p))

      if (annsum_potential_gpp(p) > 0.0) then
         avail_retransn(p) = (annmax_retransn(p)/2.0)*(gpp(p)/annsum_potential_gpp(p))/dt
      else
         avail_retransn(p) = 0.0
      end if

      ! make sure available retrans N doesn't exceed storage
      avail_retransn(p) = min(avail_retransn(p), retransn(p)/dt)

      ! modify plant N demand according to the availability of
      ! retranslocated N
      ! take from retransn pool at most the flux required to meet
      ! plant ndemand

      if (plant_ndemand(p) > avail_retransn(p)) then
         retransn_to_npool(p) = avail_retransn(p)
      else
         retransn_to_npool(p) = plant_ndemand(p)
      end if
      plant_ndemand(p) = plant_ndemand(p) - retransn_to_npool(p)

   end do ! end pft loop

   ! now use the p2c routine to get the column-averaged plant_ndemand
   allocate(col_plant_ndemand(lbc:ubc))
   call p2c(num_soilc,filter_soilc,plant_ndemand,col_plant_ndemand)

   ! column loop to resolve plant/heterotroph competition for mineral N
   do fc=1,num_soilc
      c = filter_soilc(fc)

      sum_ndemand = col_plant_ndemand(c) + potential_immob(c)

      if (sum_ndemand*dt < sminn(c)) then
         ! N availability is not limiting immobilization of plant
         ! uptake, and both can proceed at their potential rates

         nlimit = 0
         fpi(c) = 1.0
         actual_immob(c) = potential_immob(c)
         sminn_to_plant(c) = col_plant_ndemand(c)

         ! under conditions of excess N, some proportion is asusmed to
         ! be lost to denitrification, in addition to the constant
         ! proportion lost in the decomposition pathways

         sminn_to_denit_excess(c) = bdnr*((sminn(c)/dt) - sum_ndemand)
      else

         ! N availability can not satisfy the sum of immobilization and
         ! plant growth demands, so these two demands compete for available
         ! soil mineral N resource.

         nlimit = 1
         if (sum_ndemand > 0.0) then
            actual_immob(c) = (sminn(c)/dt)*(potential_immob(c) / sum_ndemand)
         else
            actual_immob(c) = 0.0
         end if

         if (potential_immob(c) > 0.0) then
            fpi(c) = actual_immob(c) / potential_immob(c)
         else
            fpi(c) = 0.0
         end if

         sminn_to_plant(c) = (sminn(c)/dt) - actual_immob(c)
      end if

      ! calculate the fraction of potential growth that can be
      ! acheived with the N available to plants

      if (col_plant_ndemand(c) > 0.0) then
         fpg(c) = sminn_to_plant(c) / col_plant_ndemand(c)
      else
         fpg(c) = 1.0
      end if

   end do ! end of column loop

   ! start new pft loop to distribute the available N between the
   ! competing pfts on the basis of relative demand, and allocate C and N to
   ! new growth and storage

   do fp=1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)

      ! set some local allocation variables
      f1 = froot_leaf(ivt(p))
      f2 = croot_stem(ivt(p))

      ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
      ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
      ! There was an error in this formula in previous version, where the coefficient
      ! was 0.004 instead of 0.0025.
      ! This variable allocation is only for trees. Shrubs have a constant
      ! allocation as specified in the pft-physiology file.  The value is also used
      ! as a trigger here: -1.0 means to use the dynamic allocation (trees).
      if (stem_leaf(ivt(p)) == -1.) then
        f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
      else
        f3 = stem_leaf(ivt(p))
      end if
      
      f4 = flivewd(ivt(p))
      g1 = grperc
      g2 = grpnow
      cnl = leafcn(ivt(p))
      cnfr = frootcn(ivt(p))
      cnlw = livewdcn(ivt(p))
      cndw = deadwdcn(ivt(p))
      fcur = fcur2(ivt(p))

      ! increase fcur linearly with ndays_active, until fcur reaches 1.0 at
      ! ndays_active = 365.  This prevents the continued storage of C and N.
      ! turning off this correction (PET, 12/11/03), instead using bgtr in
      ! phenology algorithm.
      !fcur = fcur + (1. - fcur)*lgsf(p)

      sminn_to_npool(p) = plant_ndemand(p) * fpg(c)
      plant_nalloc(p) = sminn_to_npool(p) + retransn_to_npool(p)

      ! calculate the associated carbon allocation, and the excess
      ! carbon flux that must be accounted for through downregulation

      plant_calloc(p) = plant_nalloc(p) * (c_allometry(p)/n_allometry(p))
      excess_cflux(p) = availc(p) - plant_calloc(p)

      ! reduce gpp fluxes due to N limitation
      if (gpp(p) > 0.0) then
         downreg(p) = excess_cflux(p)/gpp(p)
         psnsun_to_cpool(p) = psnsun_to_cpool(p)*(1. - downreg(p))
         psnshade_to_cpool(p) = psnshade_to_cpool(p)*(1. - downreg(p))
      end if

      ! calculate the amount of new leaf C dictated by these allocation
      ! decisions, and calculate the daily fluxes of C and N to current
      ! growth and storage pools

      ! fcur is the proportion of this day's growth that is displayed now,
      ! the remainder going into storage for display next year through the
      ! transfer pools

      nlc = plant_calloc(p) / c_allometry(p)
      cpool_to_leafc(p)          = nlc * fcur
      cpool_to_leafc_storage(p)  = nlc * (1. - fcur)
      cpool_to_frootc(p)         = nlc * f1 * fcur
      cpool_to_frootc_storage(p) = nlc * f1 * (1. - fcur)
      if (woody(ivt(p)) == 1.) then
         cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
         cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1. - fcur)
         cpool_to_deadstemc(p)          = nlc * f3 * (1. - f4) * fcur
         cpool_to_deadstemc_storage(p)  = nlc * f3 * (1. - f4) * (1. - fcur)
         cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
         cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1. - fcur)
         cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1. - f4) * fcur
         cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1. - f4) * (1. - fcur)
      end if

      ! corresponding N fluxes
      npool_to_leafn(p)          = (nlc / cnl) * fcur
      npool_to_leafn_storage(p)  = (nlc / cnl) * (1. - fcur)
      npool_to_frootn(p)         = (nlc * f1 / cnfr) * fcur
      npool_to_frootn_storage(p) = (nlc * f1 / cnfr) * (1. - fcur)
      if (woody(ivt(p)) == 1.) then
         npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw) * fcur
         npool_to_livestemn_storage(p)  = (nlc * f3 * f4 / cnlw) * (1. - fcur)
         npool_to_deadstemn(p)          = (nlc * f3 * (1. - f4) / cndw) * fcur
         npool_to_deadstemn_storage(p)  = (nlc * f3 * (1. - f4) / cndw) * (1. - fcur)
         npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
         npool_to_livecrootn_storage(p) = (nlc * f2 * f3 * f4 / cnlw) * (1. - fcur)
         npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1. - f4) / cndw) * fcur
         npool_to_deadcrootn_storage(p) = (nlc * f2 * f3 * (1. - f4) / cndw) * (1. - fcur)
      end if

      ! Calculate the amount of carbon that needs to go into growth
      ! respiration storage to satisfy all of the storage growth demands.
      ! Allows for the fraction of growth respiration that is released at the
      ! time of fixation, versus the remaining fraction that is stored for
      ! release at the time of display. Note that all the growth respiration
      ! fluxes that get released on a given timestep are calculated in growth_resp(),
      ! but that the storage of C for growth resp during display of transferred
      ! growth is assigned here.

      gresp_storage = cpool_to_leafc_storage(p) + cpool_to_frootc_storage(p)
      if (woody(ivt(p)) == 1.) then
         gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
         gresp_storage = gresp_storage + cpool_to_deadstemc_storage(p)
         gresp_storage = gresp_storage + cpool_to_livecrootc_storage(p)
         gresp_storage = gresp_storage + cpool_to_deadcrootc_storage(p)
      end if
      cpool_to_gresp_storage(p) = gresp_storage * g1 * (1. - g2)

   end do ! end pft loop

   deallocate(col_plant_ndemand)

end subroutine CNAllocation

end module CNAllocationMod
