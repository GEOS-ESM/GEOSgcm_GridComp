module CNEcosystemDynMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNEcosystemDynMod
!
! !DESCRIPTION:
! Ecosystem dynamics: phenology, vegetation
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNEcosystemDyn   ! Ecosystem dynamics: phenology, vegetation
!
! !REVISION HISTORY:
! Created by Peter Thornton
! 19 May 2009: PET - modified to include call to harvest routine
!
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !PRIVATE TYPES:
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNEcosystemDyn
!
! !INTERFACE:
  subroutine CNEcosystemDyn(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                     num_soilp, filter_soilp, doalb, spin)
!
! !DESCRIPTION:
! The core CN code is executed here. Calculates fluxes for maintenance
! respiration, decomposition, allocation, phenology, and growth respiration.
! These routines happen on the radiation time step so that canopy structure
! stays synchronized with albedo calculations.
!
! !USES:
    use clmtype
!   use spmdMod              , only: masterproc
    use CNSetValueMod        , only: CNZeroFluxes
    use CNNDynamicsMod       , only: CNNDeposition,CNNFixation, CNNLeaching
    use CNMRespMod           , only: CNMResp
    use CNDecompMod          , only: CNDecompAlloc
    use CNPhenologyMod       , only: CNPhenology
    use CNGRespMod           , only: CNGResp
    use CNCStateUpdate1Mod   , only: CStateUpdate1,CStateUpdate0
    use CNNStateUpdate1Mod   , only: NStateUpdate1
    use CNGapMortalityMod    , only: CNGapMortality
    use CNCStateUpdate2Mod   , only: CStateUpdate2, CStateUpdate2h
    use CNNStateUpdate2Mod   , only: NStateUpdate2, NStateUpdate2h
    use CNFireMod            , only: CNFireArea, CNFireFluxes
    use CNCStateUpdate3Mod   , only: CStateUpdate3
    use CNNStateUpdate3Mod   , only: NStateUpdate3
    use CNPrecisionControlMod, only: CNPrecisionControl
    use CNVegStructUpdateMod , only: CNVegStructUpdate
!   use CNAnnualUpdateMod    , only: CNAnnualUpdate
    use CNSummaryMod         , only: CSummary, NSummary
    use CNHarvestMod         , only: CNHarvest
    use CNWoodProductsMod    , only: CNWoodProducts
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc        ! column bounds
    integer, intent(in) :: lbp, ubp        ! pft bounds
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(ubp-lbp+1) ! filter for soil pfts
    logical, intent(in) :: doalb           ! true = surface albedo calculation time step
    logical, intent(in) :: spin            ! true = spinup for vegetation
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 10/22/03, Peter Thornton: created from EcosystemDyn during migration to
!                           new vector code.
! 11/3/03, Peter Thornton: removed update of elai, esai, frac_veg_nosno_alb.
!     These are now done in CNVegStructUpdate(), which is called
!     prior to SurfaceAlbedo().
! 11/13/03, Peter Thornton: switched from nolake to soil filtering.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
! local pointers to implicit out arguments
!
! !OTHER LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

       ! Call the main CN routines
       call CNZeroFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNNDeposition(lbc, ubc)

       call CNNFixation(num_soilc,filter_soilc)

       call CNMResp(lbc, ubc, num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNDecompAlloc(lbp, ubp, lbc, ubc, num_soilc, filter_soilc, num_soilp, filter_soilp, spin)

       ! CNphenology needs to be called after CNdecompAlloc, becuase it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call CNPhenology(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNGResp(num_soilp, filter_soilp)
       
       call CStateUpdate0(num_soilp, filter_soilp)

       call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNGapMortality(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
!!!    call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       call CNWoodProducts(num_soilc, filter_soilc)
       
       call CNFireArea(num_soilc, filter_soilc)

       call CNFireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNNLeaching(lbc, ubc, num_soilc, filter_soilc)

       call CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)

    if (doalb) then   
       call CNVegStructUpdate(num_soilp, filter_soilp)
    end if

!       call CNAnnualUpdate(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       call CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       call NSummary(num_soilc, filter_soilc, num_soilp, filter_soilp)

  end subroutine CNEcosystemDyn
!-----------------------------------------------------------------------
end  module CNEcosystemDynMod
