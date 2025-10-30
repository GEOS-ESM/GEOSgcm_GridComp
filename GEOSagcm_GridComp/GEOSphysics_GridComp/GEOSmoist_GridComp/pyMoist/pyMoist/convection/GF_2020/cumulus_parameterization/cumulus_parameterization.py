from ndsl import StencilFactory, QuantityFactory, ndsl_log
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
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
from pyMoist.convection.GF_2020.cumulus_parameterization.sounding.sounding import Sounding, GATESounding
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
from pyMoist.convection.GF_2020.cumulus_parameterization.parcel_moist_static_energy.parcel_moist_static_energy import (
    ParcelMoistStaticEnergy,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.entrainment.entrainment import (
    VerticalEntrainmentRate,
    DowndraftEntrainmentProfiles,
    StableDetrainment,
    CalculateMassEntrainmentDetrainment,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.moist_static_energy.moist_static_energy import (
    UpdateMoistStaticEnergy,
    FirstGuessMoistStaticEnergy,
    MoistStaticEnergyInsideCloud,
)

from pyMoist.convection.GF_2020.cumulus_parameterization.buoyancy.buoyancy import GetBuoyancy
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
from pyMoist.convection.GF_2020.cumulus_parameterization.diurnal_cycle.diurnal_cycle import DiurnalCycle
from pyMoist.convection.GF_2020.cumulus_parameterization.cape_removal.cape_removal import CAPERemoval
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
from pyMoist.convection.GF_2020.cumulus_parameterization.feedback.feedback import Feedback
from pyMoist.convection.GF_2020.cumulus_parameterization.prepare_output.prepare_output import PrepareOutput


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

        self._get_lcl = GetLCL()

        self._parcel_moist_static_energy = ParcelMoistStaticEnergy()

        self.vertical_entrainment_rate = VerticalEntrainmentRate()

        self._convective_cloud_base_level = ConvectiveCloudBaseLevel()

        self._downdraft_entraiment_profiles = DowndraftEntrainmentProfiles()

        self._update_moist_static_energy = UpdateMoistStaticEnergy()

        self._stable_detrainment = StableDetrainment()

        self._cloud_top = CloudTop()

        self._updraft_mass_flux_profile = UpdraftMassFluxProfile()

        self._calculate_mass_entrainment_detrainment = CalculateMassEntrainmentDetrainment()

        self._first_guess_moist_static_energy = FirstGuessMoistStaticEnergy()

        self._get_buoyancy = GetBuoyancy()

        self._c1d_profile = C1DProfile()

        self._updraft_moisture_profile = UpdraftMoistureProfile()

        self._melting_profile = MeltingProfile()

        self._moist_static_energy_and_momentum_budget = UpdraftMoistStaticEnergyAndMomentumBudget()

        self._in_cloud_updraft_air_temperature = UpdraftInCloudUpdraftAirTemperature()

        self._vertical_velosity = VerticalVelosity()

        self._downdraft_origin_level = DowndraftOriginLevel()

        self._downdraft_normalized_mass_flux = DowndraftNormalizedMassFlux()

        self._downdraft_lateral_mass_flux = DowndraftLateralMassFlux()

        self._downdraft_wet_bulb = DowndraftWetBlub()

        self._downdraft_moist_static_energy_and_moisture_budget = (
            DowndraftMoistStaticEnergyAndMoistureBudget()
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
                self._environment_conditions(
                    state=state,
                    locals=locals,
                    saturation_tables=saturation_tables,
                    plume_dependent_constants=self.plume_dependent_constants,
                    data_type=0,
                )
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
                self._environment_cloud_levels(
                    state=state,
                    locals=locals,
                    saturation_tables=saturation_tables,
                    plume_dependent_constants=self.plume_dependent_constants,
                    data_type=0,
                )
                self._environment_cloud_levels(
                    state=state,
                    locals=locals,
                    saturation_tables=saturation_tables,
                    plume_dependent_constants=self.plume_dependent_constants,
                    data_type=1,
                )

                # get air density at full layer (model levels) by hydrostatic balance (kg/m3)
                self._hydrostatic_air_density(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # partition between liq/ice cloud contents
                self._partition_liquid_ice(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # get the maximum origin level up an updraft
                self._maximum_updraft_origin_level(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # level where detrainment for downdraft starts
                self._downdraft_detrainment_level(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # determine level with highest moist static energy content (k_max_mse)
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
                self._get_lcl()

                # determine the moist static energy of air parcels at source level
                self._parcel_moist_static_energy()

                # determine the vertical entrainment/detrainment rates
                self._environment_conditions()

                # determine level of convective cloud base
                self._convective_cloud_base_level()

                # define entrainment/detrainment profiles for downdrafts
                self._downdraft_entraiment_profiles()

                # update unforced & forced moist static energy
                self._update_moist_static_energy()

                # increase detrainment in stable layers
                self._stable_detrainment()

                # use cloud for plumes
                self._cloud_top()

                # determine the normalized mass flux profile for updraft
                self._updraft_mass_flux_profile()

                # calculate mass entrainment and detrainment
                self._calculate_mass_entrainment_detrainment()

                # 1st guess for moist static energy
                self._first_guess_moist_static_energy()

                # Get buoyancy of updrafts
                self._get_buoyancy()

                # get "c1d" profile
                self._c1d_profile()

                # get updraft profile
                self._updraft_moisture_profile()

                # get melting profile
                self._melting_profile()

                # updraft moist static energy + momentum budget
                self._moist_static_energy_and_momentum_budget()

                # Get buoyancy of updrafts
                self._get_buoyancy()
                self._get_buoyancy()

                # calculate in-cloud/updraft air temperature for vertical velocity
                self._in_cloud_updraft_air_temperature()

                # vertical velocity
                self._vertical_velosity()

                # downdraft origin level
                self._downdraft_origin_level()

                # downdraft normalized mass flux
                self._downdraft_normalized_mass_flux()

                # lateral mass fluxes associated with downdrafts
                self._downdraft_lateral_mass_flux()

                # wet bulb temperature and moisture at downdraft origin level
                self._downdraft_wet_bulb()

                # downdraft moist static energy + moisture budget
                self._downdraft_moist_static_energy_and_moisture_budget()

                # calculate moisture properties of downdraft
                self._downdraft_moisture_properties()

                # calculate workfunctions for updrafts
                self._updraft_initial_workfunctions()

                # calculate CIN for updrafts
                self._updraft_cin()

                # trigger function: KE+CIN < 0 --> no convection
                self._trigger_function_convection()

                # calculate in-cloud/updraft and downdraft air temperature for vertical velocity
                self._in_cloud_temperature()

                # diurnal cycle section
                self._diurnal_cycle()

                # Bechtold et al 2008 time-scale of cape removal
                self._cape_removal()

                # Trigger function based on Xie et al 2019
                self._trigger_function_xie()

                # determine downdraft strength in terms of windshear
                self._downdraft_windshear()

                # get the environmental mass flux
                self._environment_mass_flux()

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
                self._moist_static_energy_inside_cloud()

                # workfunctions for updraft
                self._updraft_update_workfunctions()

                # large scale forcing
                # calculate cloud base mass flux
                self._cloud_base_mass_flux()

                # Include kinetic energy dissipation converted to heating
                self._kinetic_energy_to_heating()

                # feedback
                self._feedback()

                # net precipitation flux (after downdraft evaporation)
                self._precipitation_flux()

                # rainfall evap below cloud base
                self._rain_evap_below_cloud_base()

                # includes effects of the remained cloud dissipation into the enviroment
                self._cloud_dissapation()

                # total (deep+mid) evaporation flux for output (units kg/kg/s)
                self._output_evaporation_flux()

                # lightning flashes density (parameterization from Lopez 2016, MWR)
                self._lightning_flash_density()

                # output precipitation (only deep plume)
                self._output_deep_precipitation()

                # for tracer convective transport / outputs
                self._tracer_output()

                # convert mass fluxes, etc...
                self._prepare_output()

                # outputs a model sounding for the stand-alone code (part 2)
                self._sounding()

                # section for atmospheric composition
                self._atmospheric_composition()

                # begin: for GATE soundings
                # NOTE not needed right now, probably will never be implemented
                self._gate_sounding()
