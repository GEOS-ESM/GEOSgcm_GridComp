from ndsl import StencilFactory, QuantityFactory, ndsl_log
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import (
    GF2020CumulusParameterizationConfig,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.state import (
    GF2020CumulusParameterizationState,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import (
    GF2020CumulusParameterizationLocals,
)
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.setup.setup import Setup
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.environment.environment import (
    EnvironmentConditions,
    EnvironmentCloudLevels,
    EnvironmentMassFlux,
    EnvironmentalSubsidence,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.sounding.sounding import (
    Sounding,
    GATESounding,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.air_density.air_density import (
    HydrostaticAirDensity,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.precip.precip import (
    PartitionLiquidIce,
    PrecipFactor,
    PrecipitationFlux,
    RainEvapBelowCloudBase,
    CloudDissipation,
    OutputEvaporationFlux,
    LightningFlassDensity,
    OutputDeepPrecipitation,
    UpdateWorkfunctionsAndCondensates,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels.get_levels import (
    MaximumUpdraftOriginLevel,
    DowndraftDetrainmentLevel,
    HighestMoistStaticEnergyLevel,
    ConvectiveCloudBaseLevel,
    GetLCL,
    CloudTop,
    DowndraftOriginLevel,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.convective_tracers.convective_tracers import (
    ColdPoolParameterization,
    TracerOutput,
    AtmosphericComposition,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.entrainment.entrainment import (
    EntrainmentRates,
    DowndraftEntrainmentProfiles,
    StableDetrainment,
    CalculateMassEntrainmentDetrainment,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.moist_static_energy.moist_static_energy import (
    ParcelMoistStaticEnergy,
    UpdateMoistStaticEnergy,
    FirstGuessMoistStaticEnergy,
    MoistStaticEnergyInsideCloud,
)

from pyMoist.convection.GF_2020.cumulus_parameterization.buoyancy.buoyancy import (
    GetBuoyancy,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.profiles.profiles import (
    C1DProfile,
    MeltingProfile,
    InCloudTemperature,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.updraft.updraft import (
    UpdraftMassFluxProfile,
    UpdraftMoistureProfile,
    UpdraftMoistStaticEnergyAndMomentumBudget,
    UpdraftInCloudUpdraftAirTemperature,
    UpdraftInitialWorkfunctions,
    UpdraftCIN,
    UpdraftUpdateWorkfunctions,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.vertical_velosity.vertical_velosity import (
    VerticalVelosity,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.downdraft.downdraft import (
    DowndraftNormalizedMassFlux,
    DowndraftLateralMassFlux,
    DowndraftWetBlub,
    DowndraftMoistStaticEnergyAndMoistureBudget,
    DowndraftMoistureProperties,
    DowndraftWindshear,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.trigger_function.trigger_function import (
    TriggerFunctionConvection,
    TriggerFunctionXie,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.diurnal_cycle.diurnal_cycle import (
    DiurnalCycle,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.cape_removal.cape_removal import (
    CAPERemoval,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.mass_conservation.mass_conservation import (
    MassConservation,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.vertical_discretization.vertical_discretization import (
    VerticalDiscretization,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.smoothing.smoothing import (
    MakeSmoother,
    ApplySmoother,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.cloud_base_mass_flux.cloud_base_mass_flux import (
    CloudBaseMassFlux,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.kinetic_energy_to_heating.kinetic_energy_to_heating import (
    KineticEnergyToHeating,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.feedback.feedback import (
    Feedback,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.prepare_output.prepare_output import (
    PrepareOutput,
)


class CumulusParameterization:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # get config from the rest of the model
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # initalize the plume dependent constants, to be set for each plume within the loop (in Setup)
        self.plume_dependent_constants = GF2020PlumeDependentConstants()

        # initalize all the subclasses
        self._setup = Setup(
            stencil_factory,
            quantity_factory,
            config,
            cumulus_parameterization_config,
        )

        self._environment_conditions = EnvironmentConditions(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._sounding = Sounding()

        self._environment_cloud_levels = EnvironmentCloudLevels(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._hydrostatic_air_density = HydrostaticAirDensity(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._partition_liquid_ice = PartitionLiquidIce(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._maximum_updraft_origin_level = MaximumUpdraftOriginLevel(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._downdraft_detrainment_level = DowndraftDetrainmentLevel(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._highest_moist_static_energy_level = HighestMoistStaticEnergyLevel(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._precip_factor = PrecipFactor()

        self._cold_pool_parameterization = ColdPoolParameterization()

        self._get_lcl = GetLCL(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._parcel_moist_static_energy = ParcelMoistStaticEnergy(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._entrainment_rates = EntrainmentRates(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._convective_cloud_base_level = ConvectiveCloudBaseLevel(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._downdraft_entraiment_profiles = DowndraftEntrainmentProfiles(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._stable_detrainment = StableDetrainment(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._cloud_top = CloudTop(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._updraft_mass_flux_profile = UpdraftMassFluxProfile(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._calculate_mass_entrainment_detrainment = (
            CalculateMassEntrainmentDetrainment(
                stencil_factory=stencil_factory,
                quantity_factory=quantity_factory,
                config=config,
                cumulus_parameterization_config=cumulus_parameterization_config,
            )
        )

        self._first_guess_moist_static_energy = FirstGuessMoistStaticEnergy(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._get_buoyancy = GetBuoyancy(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._c1d_profile = C1DProfile()

        self._updraft_moisture_profile = UpdraftMoistureProfile()

        self._melting_profile = MeltingProfile()

        self._moist_static_energy_and_momentum_budget = (
            UpdraftMoistStaticEnergyAndMomentumBudget()
        )

        self._in_cloud_updraft_air_temperature = UpdraftInCloudUpdraftAirTemperature(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._vertical_velosity = VerticalVelosity()

        self._downdraft_origin_level = DowndraftOriginLevel()

        self._downdraft_normalized_mass_flux = DowndraftNormalizedMassFlux()

        self._downdraft_lateral_mass_flux = DowndraftLateralMassFlux()

        self._downdraft_wet_bulb = DowndraftWetBlub()

        self._downdraft_moist_static_energy_and_moisture_budget = (
            DowndraftMoistStaticEnergyAndMoistureBudget(
                stencil_factory=stencil_factory,
                quantity_factory=quantity_factory,
                config=config,
                cumulus_parameterization_config=cumulus_parameterization_config,
            )
        )

        self._downdraft_moisture_properties = DowndraftMoistureProperties()

        self._updraft_initial_workfunctions = UpdraftInitialWorkfunctions()

        self._updraft_cin = UpdraftCIN()

        self._trigger_function_convection = TriggerFunctionConvection()

        self._in_cloud_temperature = InCloudTemperature()

        self._diurnal_cycle = DiurnalCycle()

        self._cape_removal = CAPERemoval()

        self._trigger_function_xie = TriggerFunctionXie()

        self._downdraft_windshear = DowndraftWindshear()

        self._environment_mass_flux = EnvironmentMassFlux()

        self._mass_conservation = MassConservation()

        self._vertical_discretization = VerticalDiscretization()

        self._environmental_subsidence = EnvironmentalSubsidence()

        self._make_smoother = MakeSmoother()

        self._apply_smoother = ApplySmoother()

        self._moist_static_energy_inside_cloud = MoistStaticEnergyInsideCloud()

        self._updraft_update_workfunctions = UpdraftUpdateWorkfunctions()

        self._cloud_base_mass_flux = CloudBaseMassFlux()

        self._kinetic_energy_to_heating = KineticEnergyToHeating()

        self._feedback = Feedback()

        self._precipitation_flux = PrecipitationFlux()

        self._rain_evap_below_cloud_base = RainEvapBelowCloudBase()

        self._cloud_dissapation = CloudDissipation()

        self._output_evaporation_flux = OutputEvaporationFlux()

        self._lightning_flash_density = LightningFlassDensity()

        self._output_deep_precipitation = OutputDeepPrecipitation()

        self._tracer_output = TracerOutput()

        self._prepare_output = PrepareOutput()

        self._update_workfunctions_and_condensates = UpdateWorkfunctionsAndCondensates()

        self._atmospheric_composition = AtmosphericComposition()

        self._gate_sounding = GATESounding()

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        saturation_tables: SaturationVaporPressureTable,
    ):
        if self.cumulus_parameterization_config.PLUME_ORDER == 0:
            plume_types = ["shallow", "mid", "deep"]
            # ndsl_log.log(msg="running cumulus parameterization with order shallow, mid, deep")
        elif self.cumulus_parameterization_config.PLUME_ORDER == 1:
            plume_types = ["shallow", "deep", "mid"]
            # ndsl_log.log(msg="running cumulus parameterization with order shallow, deep, mid")
        else:
            raise NotImplementedError("plume order not impelemented")

        for plume in plume_types:
            # setup constants for the current plume, reset necessary fields, prefill necessary fields
            # NOTE test F2020_CumulusParameterization_Setup_{plume}:
            # NOTE      deep ✅
            # NOTE      mid ✅
            # NOTE      shallow ❌ (python vs fortran initalization difference)
            self._setup(
                state,
                locals,
                saturation_tables,
                self.plume_dependent_constants,
                plume,
            )

            if self.plume_dependent_constants.ENABLE_PLUME == 1:
                # environmental conditions, first heights
                # calculate moist static energy, heights, environmental saturation mixing ratio
                # NOTE test GF2020_CumulusParameterization_EnvironmentConditions_1_{plume}:
                # NOTE      deep ❌ (two vars, each one point off by two ulp)
                # NOTE      mid ❌ (two vars, each one point off by two ulp)
                # NOTE      shallow ✅
                self._environment_conditions(
                    state=state,
                    locals=locals,
                    saturation_tables=saturation_tables,
                    plume_dependent_constants=self.plume_dependent_constants,
                    data_type=0,
                )
                # NOTE test GF2020_CumulusParameterization_EnvironmentConditions_2_{plume}:
                # NOTE      deep ❌ (two vars, each one point off by two ulp)
                # NOTE      mid ❌ (two vars, each one point off by two ulp)
                # NOTE      shallow ✅
                self._environment_conditions(
                    state=state,
                    locals=locals,
                    saturation_tables=saturation_tables,
                    plume_dependent_constants=self.plume_dependent_constants,
                    data_type=1,
                )

                # outputs a model sounding for the stand-alone code (part 1)
                self._sounding()
                if self.cumulus_parameterization_config.OUTPUT_SOUNDING != 0:
                    ndsl_log.error(" GF2020 output sounding not supported")

                # environmental values on cloud levels
                # NOTE test GF2020_CumulusParameterization_EnvironmentCloudLevels_1_{plume}:
                # NOTE.     deep ❌ (worst fail rate 0.01%, worse fail 32 ULP)
                # NOTE.     mid ❌ (worst fail rate 0.01%, worse fail 32 ULP)
                # NOTE.     shallow ✅
                self._environment_cloud_levels(
                    state=state,
                    locals=locals,
                    saturation_tables=saturation_tables,
                    plume_dependent_constants=self.plume_dependent_constants,
                    data_type=0,
                )
                # NOTE test GF2020_CumulusParameterization_EnvironmentCloudLevels_2_{plume}:
                # NOTE      deep ❌ (worst fail rate 0.02%, worse fail 84 ULP)
                # NOTE      mid ❌ (worst fail rate 0.02%, worse fail 84 ULP)
                # NOTE      shallow ✅
                self._environment_cloud_levels(
                    state=state,
                    locals=locals,
                    saturation_tables=saturation_tables,
                    plume_dependent_constants=self.plume_dependent_constants,
                    data_type=1,
                )

                # get air density at full layer (model levels) by hydrostatic balance (kg/m3)
                # NOTE test GF2020_CumulusParameterization_HydrostaticAirDensity_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._hydrostatic_air_density(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # partition between liq/ice cloud contents
                # NOTE test GF2020_CumulusParameterization_PartitionLiquidIceAndGetLevels_{plume}:
                # NOTE      deep ⚠️⚠️⚠️ NEEDS ATTENTION
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._partition_liquid_ice(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # get the maximum origin level up an updraft
                # NOTE test GF2020_CumulusParameterization_PartitionLiquidIceAndGetLevels_{plume}:
                # NOTE      deep ⚠️⚠️⚠️ NEEDS ATTENTION
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._maximum_updraft_origin_level(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # level where detrainment for downdraft starts
                # NOTE test GF2020_CumulusParameterization_PartitionLiquidIceAndGetLevels_{plume}:
                # NOTE      deep ⚠️⚠️⚠️ NEEDS ATTENTION
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._downdraft_detrainment_level(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # determine level with highest moist static energy content (k_max_mse)
                # NOTE test GF2020_CumulusParameterization_HighestMoistStaticEnergyLevel_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._highest_moist_static_energy_level(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # get the pickup of ensemble ave prec, following Neelin et al 2009.
                # NOTE this runs in fortran but the output is never used - so it is not implemented
                self._precip_factor()

                # cold pool parameterization and convective memory
                # NOTE CONVECTION TRACER BLOCK
                # NOTE not called in experiment used to design this, so it is not implemented
                self._cold_pool_parameterization()

                # determine LCL for the air parcels with highest moist static energy
                # NOTE test GF2020_CumulusParameterization_GetLCL_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._get_lcl(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # determine the moist static energy of air parcels at source level
                # NOTE test GF2020_CumulusParameterization_ParcelMoistStaticEnergy_1_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._parcel_moist_static_energy(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # determine the vertical entrainment/detrainment rates
                # NOTE test GF2020_CumulusParameterization_EntrainmentRates_{plume}:
                # NOTE      deep ⚠️⚠️⚠️ NEEDS ATTENTION (maybe bad constant?)
                # NOTE      mid ⚠️⚠️⚠️ NEEDS ATTENTION (maybe bad constant?)
                # NOTE      shallow ✅
                self._entrainment_rates(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # determine level of convective cloud base
                # NOTE test GF2020_CumulusParameterization_ConvectiveCloudBaseLevel_{plume}:
                # NOTE      deep ⚠️⚠️⚠️ NEEDS ATTENTION (gtir error)
                # NOTE      mid ⚠️⚠️⚠️ TEST DOES NOT YET EXIST
                # NOTE      shallow ⚠️⚠️⚠️ TEST DOES NOT YET EXIST
                self._convective_cloud_base_level(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # define entrainment/detrainment profiles for downdrafts
                # NOTE test GF2020_CumulusParameterization_DowndraftEntrainmentProfiles_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._downdraft_entraiment_profiles(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # update unforced & forced moist static energy
                # NOTE test GF2020_CumulusParameterization_ParcelMoistStaticEnergy_2_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._parcel_moist_static_energy(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # increase detrainment in stable layers
                # NOTE test GF2020_CumulusParameterization_StableDetrainment_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._stable_detrainment(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # use cloud for plumes
                # NOTE test GF2020_CumulusParameterization_CloudTop_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._cloud_top(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # determine the normalized mass flux profile for updraft
                # NOTE test GF2020_CumulusParameterization_UpdraftMassFluxProfile_{plume}:
                # NOTE      deep ⚠️⚠️⚠️ UNFINISHED, TEST DOES NOT EXIST
                # NOTE      mid ⚠️⚠️⚠️ UNFINISHED, TEST DOES NOT EXIST
                # NOTE      shallow ⚠️⚠️⚠️ UNFINISHED, TEST DOES NOT EXIST
                self._updraft_mass_flux_profile(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # calculate mass entrainment and detrainment
                # NOTE test GF2020_CumulusParameterization_CalculateMassEntrainmentDetrainment_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._calculate_mass_entrainment_detrainment(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # 1st guess for moist static energy
                # NOTE test GF2020_CumulusParameterization_CalculateMassEntrainmentDetrainment_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ⚠️⚠️⚠️ NEEDS ATTENTION
                # NOTE      shallow ⚠️⚠️⚠️ NEEDS ATTENTION
                self._first_guess_moist_static_energy(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # Get buoyancy of updrafts
                # NOTE test GF2020_CumulusParameterization_GetBuoyancy_1_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._get_buoyancy(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # get "c1d" profile
                self._c1d_profile()

                # get updraft profile
                self._updraft_moisture_profile()

                # get melting profile
                # NOTE test GF2020_CumulusParameterization_MeltingProfile_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._melting_profile(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # updraft moist static energy + momentum budget
                self._moist_static_energy_and_momentum_budget()

                # Get buoyancy of updrafts
                # NOTE test GF2020_CumulusParameterization_GetBuoyancy_2_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._get_buoyancy(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )
                # NOTE test GF2020_CumulusParameterization_GetBuoyancy_3_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._get_buoyancy(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # calculate in-cloud/updraft air temperature for vertical velocity
                # NOTE ported, but untested
                self._in_cloud_updraft_air_temperature(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # vertical velocity
                self._vertical_velosity()

                # downdraft origin level
                self._downdraft_origin_level()

                # downdraft normalized mass flux
                self._downdraft_normalized_mass_flux()

                # lateral mass fluxes associated with downdrafts
                self._downdraft_lateral_mass_flux()

                # wet bulb temperature and moisture at downdraft origin level
                # NOTE ported, not tested
                self._downdraft_wet_bulb(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # downdraft moist static energy + moisture budget
                # NOTE ported, but untested
                self._downdraft_moist_static_energy_and_moisture_budget(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # calculate moisture properties of downdraft
                self._downdraft_moisture_properties()

                # calculate workfunctions for updrafts
                # NOTE ported, not tested (should pass)
                self._updraft_initial_workfunctions(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # calculate CIN for updrafts
                # NOTE ported, not tested (should pass)
                self._updraft_cin(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # trigger function: KE+CIN < 0 --> no convection
                self._trigger_function_convection(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # calculate in-cloud/updraft and downdraft air temperature for vertical velocity
                # NOTE ported, but untested
                self._in_cloud_temperature(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # diurnal cycle section
                # NOTE ported, but untested
                self._diurnal_cycle(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # Bechtold et al 2008 time-scale of cape removal
                self._cape_removal()

                # Trigger function based on Xie et al 2019
                self._trigger_function_xie()

                # determine downdraft strength in terms of windshear
                # NOTE ported, not tested
                self._downdraft_windshear(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # get the environmental mass flux
                # NOTE ported, but untested
                self._environment_mass_flux(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # check mass conservation
                self._mass_conservation()

                # change per unit mass that a model cloud would modify the environment
                self._vertical_discretization()

                # apply environmental subsidence on grid-scale ice and
                # liq water contents, and cloud fraction (Upwind scheme)
                self._environmental_subsidence()

                # make the smoothness procedure
                self._make_smoother()

                # using smoothed tendencies, calculate changed environmental profiles
                self._apply_smoother()

                # calculate moist static energy, heights, environmental saturation mixing ratio
                self._environment_conditions()

                # environmental values on cloud levels
                self._environment_cloud_levels()

                # static control
                # moist static energy inside cloud
                # NOTE ported, but untested
                self._moist_static_energy_inside_cloud(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # workfunctions for updraft
                self._updraft_update_workfunctions()

                # large scale forcing
                # calculate cloud base mass flux
                self._cloud_base_mass_flux()

                # Include kinetic energy dissipation converted to heating
                # NOTE ported, but untested
                self._kinetic_energy_to_heating(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # feedback
                self._feedback()

                # net precipitation flux (after downdraft evaporation)
                # NOTE ported, not tested
                self._precipitation_flux(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # rainfall evap below cloud base
                self._rain_evap_below_cloud_base()

                # includes effects of the remained cloud dissipation into the enviroment
                self._cloud_dissapation()

                # total (deep+mid) evaporation flux for output (units kg/kg/s)
                # NOTE ported, not tested
                self._output_evaporation_flux(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # lightning flashes density (parameterization from Lopez 2016, MWR)
                self._lightning_flash_density()

                # output precipitation (only deep plume)
                # NOTE ported, not tested
                self._output_deep_precipitation(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # for tracer convective transport / outputs
                # NOTE ported, not tested
                self._tracer_output(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # convert mass fluxes, etc...
                # NOTE ported, not tested
                self._prepare_output(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # outputs a model sounding for the stand-alone code (part 2)
                self._sounding()

                # section for atmospheric composition
                self._atmospheric_composition()

                # begin: for GATE soundings
                # NOTE not needed right now, probably will never be implemented
                self._gate_sounding()
