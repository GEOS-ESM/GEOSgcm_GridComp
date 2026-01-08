from ndsl import StencilFactory, QuantityFactory, ndsl_log
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
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
from pyMoist.convection.GF_2020.cumulus_parameterization.environment import (
    environment_conditions,
    environment_cloud_levels,
    environment_mass_flux,
    EnvironmentalSubsidence,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.sounding.sounding import (
    Sounding,
    GATESounding,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.air_density import (
    hydrostatic_air_density,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.precip import (
    partition_liquid_ice,
    PrecipFactor,
    PrecipitationFlux,
    RainEvapBelowCloudBase,
    CloudDissipation,
    OutputEvaporationFlux,
    LightningFlassDensity,
    OutputDeepPrecipitation,
    UpdateWorkfunctionsAndCondensates,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels import (
    find_maximum_updraft_origin_level,
    find_detrainmet_start_level,
    find_highest_moist_static_energy_level,
    find_lcl,
    set_start_level,
    get_convective_cloud_base_level,
    updraft_rates_pdf,
    cloud_top_checks,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.convective_tracers.convective_tracers import (
    ColdPoolParameterization,
    TracerOutput,
    AtmosphericComposition,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.entrainment import (
    entrainment_rates,
    downdraft_entrainment_profiles,
    unknown_find_level,
    compute_lateral_massflux,
    compute_uc_vc,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.moist_static_energy import (
    parcel_moist_static_energy,
    first_guess_moist_static_energy,
    MoistStaticEnergyInsideCloud,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import updraft_vertical_velocity
from pyMoist.convection.GF_2020.cumulus_parameterization.buoyancy import get_buoyancy
from pyMoist.convection.GF_2020.cumulus_parameterization.profiles import (
    C1DProfile,
    get_melting_profile,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.updraft import (
    UpdraftMassFlux,
    updraft_moisture,
    updraft_moist_static_energy_and_momentum_budget,
    updraft_temperature,
    UpdraftInitialWorkfunctions,
    UpdraftCIN,
    UpdraftUpdateWorkfunctions,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.triggers import (
    convection_trigger,
    XieTriggerFunction,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.downdraft import (
    DowndraftOriginLevel,
    DowndraftNormalizedMassFlux,
    downdraft_lateral_massflux,
    DowndraftWetBlub,
    downdraft_moist_static_energy_and_buoyancy,
    downdraft_moisture,
    downdraft_temperature,
    DowndraftWindShear,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.diurnal_cycle import DiurnalCycle
from pyMoist.convection.GF_2020.cumulus_parameterization.mass_conservation import MassConservation
from pyMoist.convection.GF_2020.cumulus_parameterization.vertical_discretization import (
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

        self._environment_conditions = stencil_factory.from_dims_halo(
            func=environment_conditions,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "SATURATION_CALCULATION_CHOICE": cumulus_parameterization_config.SATURATION_CALCULATION_CHOICE
            },
        )

        self._sounding = Sounding()

        self._environment_cloud_levels = stencil_factory.from_dims_halo(
            func=environment_cloud_levels,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"CLOUD_LEVEL_GRID": cumulus_parameterization_config.CLOUD_LEVEL_GRID},
        )

        self._hydrostatic_air_density = stencil_factory.from_dims_halo(
            func=hydrostatic_air_density,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._partition_liquid_ice = stencil_factory.from_dims_halo(
            func=partition_liquid_ice,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "MELT_GLAC": cumulus_parameterization_config.MELT_GLAC,
                "FRAC_MODIS": cumulus_parameterization_config.FRAC_MODIS,
            },
        )

        self._find_maximum_updraft_origin_level = stencil_factory.from_dims_halo(
            func=find_maximum_updraft_origin_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._find_detrainmet_start_level = stencil_factory.from_dims_halo(
            func=find_detrainmet_start_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._find_highest_moist_static_energy_level = stencil_factory.from_dims_halo(
            func=find_highest_moist_static_energy_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._precip_factor = PrecipFactor()

        self._cold_pool_parameterization = ColdPoolParameterization()

        self._find_lcl = stencil_factory.from_dims_halo(
            func=find_lcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
                "ADV_TRIGGER": config.ADV_TRIGGER,
            },
        )

        self._parcel_moist_static_energy = stencil_factory.from_dims_halo(
            func=parcel_moist_static_energy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
            },
        )

        self._entrainment_rates = stencil_factory.from_dims_halo(
            func=entrainment_rates,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF},
        )

        self._set_start_level = stencil_factory.from_dims_halo(
            func=set_start_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._get_convective_cloud_base_level = stencil_factory.from_dims_halo(
            func=get_convective_cloud_base_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "OVERSHOOT": cumulus_parameterization_config.OVERSHOOT,
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                "MOIST_TRIGGER": cumulus_parameterization_config.MOIST_TRIGGER,
                "USE_MEMORY": cumulus_parameterization_config.USE_MEMORY,
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
            },
        )

        self._downdraft_entraiment_profiles = stencil_factory.from_dims_halo(
            func=downdraft_entrainment_profiles,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"DOWNDRAFT": cumulus_parameterization_config.DOWNDRAFT},
        )

        self._unknown_find_level = stencil_factory.from_dims_halo(
            func=unknown_find_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._updraft_rates_pdf = stencil_factory.from_dims_halo(
            func=updraft_rates_pdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"OVERSHOOT": cumulus_parameterization_config.OVERSHOOT},
        )

        self._cloud_top_checks = stencil_factory.from_dims_halo(
            func=cloud_top_checks,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._updraft_mass_flux = UpdraftMassFlux(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._compute_lateral_massflux = stencil_factory.from_dims_halo(
            func=compute_lateral_massflux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._compute_uc_vc = stencil_factory.from_dims_halo(
            func=compute_uc_vc,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD
            },
        )

        self._first_guess_moist_static_energy = stencil_factory.from_dims_halo(
            func=first_guess_moist_static_energy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._get_buoyancy = stencil_factory.from_dims_halo(
            func=get_buoyancy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._c1d_profile = C1DProfile(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._updraft_moisture_profile = stencil_factory.from_dims_halo(
            func=updraft_moisture,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
                "USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES": cumulus_parameterization_config.USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES,
                "AUTOCONV": config.AUTOCONV,
                "CRITICAL_MIXING_RATIO_OVER_OCEAN": cumulus_parameterization_config.CRITICAL_MIXING_RATIO_OVER_OCEAN,
                "CRITICAL_MIXING_RATIO_OVER_LAND": cumulus_parameterization_config.CRITICAL_MIXING_RATIO_OVER_LAND,
                "FRAC_MODIS": cumulus_parameterization_config.FRAC_MODIS,
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
            },
        )

        self._melting_profile = stencil_factory.from_dims_halo(
            func=get_melting_profile,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "MELT_GLAC": cumulus_parameterization_config.MELT_GLAC,
            },
        )

        self._updraft_moist_static_energy_and_momentum_budget = stencil_factory.from_dims_halo(
            func=updraft_moist_static_energy_and_momentum_budget,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES": cumulus_parameterization_config.USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES,
                "PRESSURE_GRADIENT_CONSTANT": cumulus_parameterization_config.PRESSURE_GRADIENT_CONSTANT,
            },
        )

        self._updraft_temperature = stencil_factory.from_dims_halo(
            func=updraft_temperature,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._updraft_vertical_velosity = stencil_factory.from_dims_halo(
            func=updraft_vertical_velocity,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._downdraft_origin_level = DowndraftOriginLevel(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._downdraft_normalized_mass_flux = DowndraftNormalizedMassFlux()

        self._downdraft_lateral_mass_flux = stencil_factory.from_dims_halo(
            func=downdraft_lateral_massflux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._downdraft_wet_bulb = DowndraftWetBlub()

        self._downdraft_moist_static_energy_and_buoyancy = stencil_factory.from_dims_halo(
            func=downdraft_moist_static_energy_and_buoyancy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_WETBULB": cumulus_parameterization_config.USE_WETBULB,
                "PGCON": cumulus_parameterization_config.PRESSURE_GRADIENT_CONSTANT,
            },
        )

        self._downdraft_moisture = stencil_factory.from_dims_halo(
            func=downdraft_moisture,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_WETBULB": cumulus_parameterization_config.USE_WETBULB,
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                "EVAP_FIX": cumulus_parameterization_config.EVAP_FIX,
            },
        )

        self._updraft_initial_workfunctions = UpdraftInitialWorkfunctions(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._updraft_cin = UpdraftCIN(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._convection_trigger = stencil_factory.from_dims_halo(
            func=convection_trigger,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"DICYCLE": cumulus_parameterization_config.DIURNAL_CYCLE},
        )

        self._downdraft_temperature = stencil_factory.from_dims_halo(
            func=downdraft_temperature,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._diurnal_cycle = DiurnalCycle(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._Xie_trigger_function = XieTriggerFunction()

        self._downdraft_windshear = DowndraftWindShear(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._environment_mass_flux = stencil_factory.from_dims_halo(
            func=environment_mass_flux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._mass_conservation = MassConservation()

        self._vertical_discretization = VerticalDiscretization(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

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
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._environment_conditions(
                    p=state.input_output.p_forced,
                    p_surface=state.input_output.p_surface,
                    t=state.input_output.t_old,
                    vapor=state.input_output.vapor_old,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    moist_static_energy=locals.environment_moist_static_energy,
                    saturation_moist_static_energy=locals.environment_saturation_moist_static_energy,
                    saturation_mixing_ratio=locals.environment_saturation_mixing_ratio,
                    geopotential_height=locals.geopotential_height,
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )
                # NOTE test GF2020_CumulusParameterization_EnvironmentConditions_2_{plume}:
                # NOTE      deep ❌ (two vars, each one point off by two ulp)
                # NOTE      mid ❌ (two vars, each one point off by two ulp)
                # NOTE      shallow ✅
                self._environment_conditions(
                    p=state.input_output.p_forced,
                    p_surface=state.input_output.p_surface,
                    t=locals.t_new,
                    vapor=locals.vapor_forced,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    moist_static_energy=locals.environment_moist_static_energy_forced,
                    saturation_moist_static_energy=locals.environment_saturation_moist_static_energy_forced,
                    saturation_mixing_ratio=locals.environment_saturation_mixing_ratio_forced,
                    geopotential_height=state.input_output.geopotential_height_forced,
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # outputs a model sounding for the stand-alone code (part 1)
                if self.cumulus_parameterization_config.CLOUD_LEVEL_GRID != 1:
                    ndsl_log.warning(
                        " GF2020 cumulus parameterization initalized with unimplemented OUTPUT_SOUNDING option. "
                        "Output soundings are not currently available. Contact support if this tool is needed."
                    )

                # environmental values on cloud levels
                # NOTE test GF2020_CumulusParameterization_EnvironmentCloudLevels_1_{plume}:
                # NOTE.     deep ❌ (worst fail rate 0.01%, worse fail 32 ULP)
                # NOTE.     mid ❌ (worst fail rate 0.01%, worse fail 32 ULP)
                # NOTE.     shallow ✅
                self._environment_cloud_levels(
                    p=state.input_output.p_forced,
                    p_surface=state.input_output.p_surface,
                    p_cloud_levels=locals.p_cloud_levels,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    geopotential_height=locals.geopotential_height,
                    geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
                    t=state.input_output.t_old,
                    t_surface=state.input_output.t_surface,
                    t_cloud_levels=locals.t_cloud_levels,
                    vapor=state.input_output.vapor_old,
                    vapor_cloud_levels=locals.vapor_cloud_levels,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    u_cloud_levels=locals.u_cloud_levels,
                    v_cloud_levels=locals.v_cloud_levels,
                    environment_saturation_mixing_ratio=locals.environment_saturation_mixing_ratio,
                    environment_saturation_mixing_ratio_cloud_levels=locals.environment_saturation_mixing_ratio_cloud_levels,
                    environment_moist_static_energy=locals.environment_moist_static_energy,
                    environment_moist_static_energy_cloud_levels=locals.environment_moist_static_energy_cloud_levels,
                    environment_saturation_moist_static_energy=locals.environment_saturation_moist_static_energy,
                    environment_saturation_moist_static_energy_cloud_levels=locals.environment_saturation_moist_static_energy_cloud_levels,
                    gamma_cloud_levels=locals.gamma_cloud_levels,
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )
                # NOTE test GF2020_CumulusParameterization_EnvironmentCloudLevels_2_{plume}:
                # NOTE      deep ❌ (worst fail rate 0.02%, worse fail 84 ULP)
                # NOTE      mid ❌ (worst fail rate 0.02%, worse fail 84 ULP)
                # NOTE      shallow ✅
                p_3d = state.output.p_cloud_levels_forced.field[
                    :, :, :, self.plume_dependent_constants.PLUME_INDEX
                ]  # TODO this should go in a stencil
                self._environment_cloud_levels(
                    p=state.input_output.p_forced,
                    p_surface=state.input_output.p_surface,
                    p_cloud_levels=p_3d,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    geopotential_height=state.input_output.geopotential_height_forced,
                    geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels_forced,
                    t=locals.t_new,
                    t_surface=state.input_output.t_surface,
                    t_cloud_levels=locals.t_cloud_levels_forced,
                    vapor=locals.vapor_forced,
                    vapor_cloud_levels=locals.vapor_cloud_levels_forced,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    u_cloud_levels=locals.u_cloud_levels,
                    v_cloud_levels=locals.v_cloud_levels,
                    environment_saturation_mixing_ratio=locals.environment_saturation_mixing_ratio_forced,
                    environment_saturation_mixing_ratio_cloud_levels=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    environment_moist_static_energy=locals.environment_moist_static_energy_forced,
                    environment_moist_static_energy_cloud_levels=locals.environment_moist_static_energy_cloud_levels_forced,
                    environment_saturation_moist_static_energy=locals.environment_saturation_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    gamma_cloud_levels=locals.gamma_cloud_levels_forced,
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )
                state.output.p_cloud_levels_forced.field[
                    :, :, :, self.plume_dependent_constants.PLUME_INDEX
                ] = p_3d  # TODO this should go in a stencil

                # get air density at full layer (model levels) by hydrostatic balance (kg/m3)
                # NOTE test GF2020_CumulusParameterization_HydrostaticAirDensity_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._hydrostatic_air_density(
                    p=state.output.p_cloud_levels_forced,
                    geopotential_height=locals.geopotential_height_cloud_levels_forced,
                    error_code=state.output.error_code,
                    air_density=locals.hydrostatic_air_density,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # partition between liq/ice cloud contents
                # NOTE test GF2020_CumulusParameterization_PartitionLiquidIceAndGetLevels_{plume}:
                # NOTE      deep ❌ 1 var failing (500 ULP) -> ice_fraction function causing issues
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._partition_liquid_ice(
                    t=locals.t_new,
                    p=state.output.p_cloud_levels_forced,
                    geopotential_height=locals.geopotential_height_cloud_levels_forced,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    surface_type=state.input.surface_type,
                    convection_fraction=state.input.convection_fraction,
                    error_code=state.output.error_code,
                    melting_layer=locals.melting_layer,
                    part_liquid_ice=locals.partition_liquid_ice,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                self._find_maximum_updraft_origin_level(
                    geopotential_height=locals.geopotential_height_cloud_levels_forced,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    error_code=state.output.error_code,
                    maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
                    MAX_UPDRAFT_ORIGIN_HEIGHT=self.plume_dependent_constants.MAX_UPDRAFT_ORIGIN_HEIGHT,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                self._find_detrainmet_start_level(
                    geopotential_height=locals.geopotential_height_cloud_levels_forced,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    error_code=state.output.error_code,
                    detrainment_start_level=locals.detrainment_start_level,
                    DETRAINMENT_CRITICAL_DEPTH=self.plume_dependent_constants.DETRAINMENT_CRITICAL_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # determine level with highest moist static energy content (k_max_mse)
                # NOTE test GF2020_CumulusParameterization_HighestMoistStaticEnergyLevel_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._find_highest_moist_static_energy_level(
                    moist_static_energy=locals.environment_moist_static_energy_cloud_levels_forced,
                    error_code=state.output.error_code,
                    maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
                    updraft_origin_level=state.output.updraft_origin_level,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
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
                self._find_lcl(
                    p=state.input_output.p_forced,
                    p_cloud_levels=locals.p_cloud_levels,
                    t_excess=locals.t_excess,
                    t_cloud_levels_forced=locals.t_cloud_levels,
                    t_perturbation=state.output.t_perturbation,
                    vapor_excess=locals.vapor_excess,
                    vapor_cloud_levels_forced=locals.vapor_cloud_levels,
                    omega=state.input_output.omega,
                    air_density=state.input_output.air_density,
                    geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    ocean_fraction=state.input.ocean_fraction,
                    updraft_origin_level=state.output.updraft_origin_level,
                    grid_length=state.input_output.grid_length,
                    lcl_level=state.output.lcl_level,
                    error_code=state.output.error_code,
                    AVERAGE_LAYER_DEPTH=self.plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # determine the moist static energy of air parcels at source level
                # NOTE test GF2020_CumulusParameterization_ParcelMoistStaticEnergy_1_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._parcel_moist_static_energy(
                    error_code=state.output.error_code,
                    t_excess=locals.t_excess,
                    vapor_excess=locals.vapor_excess,
                    add_buoyancy=locals.add_buoyancy,
                    ocean_fraction=state.input.ocean_fraction,
                    updraft_origin_level=state.output.updraft_origin_level,
                    p=state.input_output.p_forced,
                    environmenet_moist_static_energy=locals.environment_moist_static_energy_cloud_levels,
                    environmenet_moist_static_energy_forced=locals.environment_moist_static_energy_cloud_levels_forced,
                    t_perturbation=state.output.t_perturbation,
                    moist_static_energy_origin_level=locals.moist_static_energy_origin_level,
                    moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
                    AVERAGE_LAYER_DEPTH=self.plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # determine the vertical entrainment/detrainment rates
                # NOTE test GF2020_CumulusParameterization_EntrainmentRates_{plume}:
                # NOTE      deep ❌ (worst fail 2.41% - max 4 ULP)
                # NOTE      mid ❌ (worst fail 0.88% - max 2 ULP)
                # NOTE      shallow ✅
                self._entrainment_rates(
                    vapor=locals.vapor_cloud_levels_forced,
                    environment_saturation_mixing_ratio=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    lcl_level=state.output.lcl_level,
                    error_code=state.output.error_code,
                    entrainment_rate=state.output.entrainment_rate,
                    detrainment_function_updraft=locals.detrainment_function_updraft,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # determine level of convective cloud base
                # NOTE test GF2020_CumulusParameterization_ConvectiveCloudBaseLevel_{plume}:
                # NOTE      deep ⚠️⚠️⚠️ MATH ERROR. major flaw in get_convective_cloud_base
                # NOTE      mid ⚠️⚠️⚠️ MATH ERROR. major flaw in get_convective_cloud_base
                # NOTE      shallow ✅
                self._set_start_level(
                    lcl_level=state.output.lcl_level,
                    start_level=locals.start_level,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                self._get_convective_cloud_base_level(
                    error_code=state.output.error_code,
                    lcl_level=state.output.lcl_level,
                    cloud_moist_static_energy_forced_transported=locals.cloud_moist_static_energy_forced_transported,
                    cap_max=locals.cap_max,
                    updraft_origin_level=state.output.updraft_origin_level,
                    start_level=locals.start_level,
                    moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
                    negative_buoyancy_depth=locals.negative_buoyancy_depth,
                    frh_lfc=locals.frh_lfc,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    entrainment_rate=state.output.entrainment_rate,
                    environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    t_excess=locals.t_excess,
                    vapor_excess=locals.vapor_excess,
                    add_buoyancy=locals.add_buoyancy,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    vapor_forced=locals.vapor_forced,
                    environment_saturation_mixing_ratio_forced=locals.environment_saturation_mixing_ratio_forced,
                    ocean_fraction=locals.ocean_fraction,
                    cap_max_increment=locals.cap_max_increment,
                    t_perturbation=state.output.t_perturbation,
                    p_forced=state.input_output.p_forced,
                    cloud_top_level=state.output.cloud_top_level,
                    AVERAGE_LAYER_DEPTH=self.plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # define entrainment/detrainment profiles for downdrafts
                # NOTE test GF2020_CumulusParameterization_DowndraftEntrainmentProfiles_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._downdraft_entraiment_profiles(
                    lateral_entrainment_rate=state.input.lateral_entrainment_rate,
                    entrainment_rate_downdraft=locals.entrainment_rate_downdraft,
                    detrainment_function_downdraft=locals.detrainment_function_downdraft,
                    scale_dependence_factor_downdraft=locals.scale_dependence_factor_downdraft,
                    plume_entrainment_rate=self.plume_dependent_constants.ENTRAINMENT_RATE,
                )

                # update unforced & forced moist static energy
                # NOTE test GF2020_CumulusParameterization_ParcelMoistStaticEnergy_2_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._parcel_moist_static_energy(
                    error_code=state.output.error_code,
                    t_excess=locals.t_excess,
                    vapor_excess=locals.vapor_excess,
                    add_buoyancy=locals.add_buoyancy,
                    ocean_fraction=state.input.ocean_fraction,
                    updraft_origin_level=state.output.updraft_origin_level,
                    p=state.input_output.p_forced,
                    environmenet_moist_static_energy=locals.environment_moist_static_energy_cloud_levels,
                    environmenet_moist_static_energy_forced=locals.environment_moist_static_energy_cloud_levels_forced,
                    t_perturbation=state.output.t_perturbation,
                    moist_static_energy_origin_level=locals.moist_static_energy_origin_level,
                    moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
                    AVERAGE_LAYER_DEPTH=self.plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # increase detrainment in stable layers
                # NOTE test GF2020_CumulusParameterization_StableDetrainment_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._unknown_find_level(
                    array=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    start_index=state.output.updraft_lfc_level,
                    end_index=locals.kstabm,
                    out_index=state.output.kstabi,
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # use cloud for plumes
                # NOTE test GF2020_CumulusParameterization_CloudTop_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._updraft_rates_pdf(
                    entrainment_rate=state.output.entrainment_rate,
                    moist_static_energy=locals.environment_moist_static_energy_forced,
                    saturation_moist_static_energy=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    moist_static_energy_origin_level=locals.moist_static_energy_origin_level_forced,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    geopotential_height=locals.geopotential_height_cloud_levels_forced,
                    cloud_moist_static_energy=locals.cloud_moist_static_energy_forced_transported,
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                self._cloud_top_checks(
                    cloud_top_level=state.output.cloud_top_level,
                    p=state.output.p_cloud_levels_forced,
                    geopotential_height=locals.geopotential_height_cloud_levels,
                    error_code=state.output.error_code,
                    last_error_code=state.input.last_error_code,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    MINIMUM_DEPTH=self.plume_dependent_constants.MINIMUM_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # determine the normalized mass flux profile for updraft
                # NOTE test GF2020_CumulusParameterization_UpdraftMassFlux_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._updraft_mass_flux(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # calculate mass entrainment and detrainment
                # NOTE test GF2020_CumulusParameterization_CalculateMassEntrainmentDetrainment_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._compute_lateral_massflux(
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    geopotential_height=locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_updraft=state.output.normalized_massflux_updraft_forced,
                    detrainment_function_updraft=locals.detrainment_function_updraft,
                    entrainment_rate=state.output.entrainment_rate,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_entrainment_updraft=locals.mass_entrainment_updraft,
                    mass_detrainment_updraft=locals.mass_detrainment_updraft,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    updraft_origin_level=state.output.updraft_origin_level,
                    pbl_level=state.input_output.pbl_level,
                    mass_entrainment_u_updraft=locals.mass_entrainment_u_updraft,
                    mass_detrainment_u_updraft=locals.mass_detrainment_u_updraft,
                    LAMBDA_DEEP=self.plume_dependent_constants.LAMBDA_DEEP,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                self._compute_uc_vc(
                    u_c=locals.u_c,
                    v_c=locals.v_c,
                    cloud_moist_static_energy=locals.cloud_moist_static_energy,
                    cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                    error_code=state.output.error_code,
                    start_level=locals.start_level,
                    moist_static_energy_origin_level=locals.moist_static_energy_origin_level,
                    moist_static_energy_origin_level_forced=locals.moist_static_energy_origin_level_forced,
                    u_cloud_levels=locals.u_cloud_levels,
                    v_cloud_levels=locals.v_cloud_levels,
                    p=state.input_output.p_forced,
                    updraft_origin_level=state.output.updraft_origin_level,
                    ocean_fraction=state.input.ocean_fraction,
                    AVERAGE_LAYER_DEPTH=self.plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # 1st guess for moist static energy
                # NOTE test GF2020_CumulusParameterization_FirstGuessMoistStaticEnergy_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._first_guess_moist_static_energy(
                    error_code=state.output.error_code,
                    start_level=locals.start_level,
                    cloud_top_level=state.output.cloud_top_level,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                    normalized_massflux_updraft=locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                    vapor_excess=locals.vapor_excess,
                    t_excess=locals.t_excess,
                    add_buoyancy=locals.add_buoyancy,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # Get buoyancy of updrafts
                # NOTE test GF2020_CumulusParameterization_GetBuoyancy_1_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._get_buoyancy(
                    lcl_level=state.output.lcl_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    cloud_moist_static_energy=locals.cloud_moist_static_energy_forced,
                    environment_moist_static_energy=locals.environment_moist_static_energy_cloud_levels_forced,
                    environment_saturation_moist_static_energy=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    buoyancy=locals.buoyancy,
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # get "c1d" profile
                # NOTE test GF2020_CumulusParameterization_C1DProfile_{plume}:
                # NOTE      deep ⚠️⚠️⚠️ DOES NOT EXECUTE IN CURRENT SIMULATION
                # NOTE      mid ⚠️⚠️⚠️ DOES NOT EXECUTE IN CURRENT SIMULATION
                # NOTE      shallow ⚠️⚠️⚠️ DOES NOT EXECUTE IN CURRENT SIMULATION
                # NOTE UNFINISHED - MANUALLY DISABLED (see class docstring)
                if False:
                    self._c1d_profile(
                        state=state,
                        locals=locals,
                        plume_dependent_constants=self.plume_dependent_constants,
                    )

                # calculate moisture properties of updraft
                # NOTE test GF2020_CumulusParameterization_UpdraftMoisture_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._updraft_moisture_profile(
                    start_level=locals.start_level,
                    error_code=state.output.error_code,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    cloud_total_water_after_entrainment_forced=locals.cloud_total_water_after_entrainment_forced,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                    cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                    updraft_column_temperature_forced=locals.updraft_column_temperature_forced,
                    ocean_fraction=state.input.ocean_fraction,
                    convection_fraction=state.input.convection_fraction,
                    surface_type=state.input.surface_type,
                    p_forced=state.input_output.p_forced,
                    cloud_top_level=state.output.cloud_top_level,
                    d_buoyancy_forced=locals.d_buoyancy_forced,
                    cloud_liquid_before_rain_forced=locals.cloud_liquid_before_rain_forced,
                    t_cloud_levels=locals.t_cloud_levels,
                    vapor_forced=locals.vapor_forced,
                    gamma_cloud_levels_forced=locals.gamma_cloud_levels_forced,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    environment_saturation_mixing_ratio_cloud_levels_forced=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    updraft_origin_level=state.output.updraft_origin_level,
                    vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
                    vapor_excess=locals.vapor_excess,
                    ccn=state.input_output.ccn,
                    mass_entrainment_updraft=locals.mass_entrainment_updraft,
                    mass_detrainment_updraft=locals.mass_detrainment_updraft,
                    psum=locals.psum,
                    psumh=locals.psumh,
                    c1d=locals.c1d,
                    add_buoyancy=locals.add_buoyancy,
                    vertical_velocity_3d=locals.vertical_velocity_3d,
                    C0=self.plume_dependent_constants.C0,
                    AVERAGE_LAYER_DEPTH=self.plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # get melting profile
                # NOTE test GF2020_CumulusParameterization_MeltingProfile_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._melting_profile(
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                    melting_layer=locals.melting_layer,
                    partition_liquid_ice=locals.partition_liquid_ice,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    melting=locals.melting,
                )

                # updraft moist static energy + momentum budget
                # NOTE test GF2020_CumulusParameterization_UpdraftMoistStaticEnergyAndMomentumBudget_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._updraft_moist_static_energy_and_momentum_budget(
                    error_code=state.output.error_code,
                    start_level=locals.start_level,
                    cloud_top_level=state.output.cloud_top_level,
                    p_forced=state.input_output.p_forced,
                    environment_moist_static_energy=locals.environment_moist_static_energy,
                    environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                    environment_moist_static_energy_cloud_levels=locals.environment_moist_static_energy_cloud_levels,
                    environment_moist_static_energy_cloud_levels_forced=locals.environment_moist_static_energy_cloud_levels_forced,
                    environment_saturation_moist_static_energy_cloud_levels=locals.environment_saturation_moist_static_energy_cloud_levels,
                    environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy=locals.cloud_moist_static_energy,
                    cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                    normalized_massflux_updraft=locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    mass_entrainment_updraft=locals.mass_entrainment_updraft,
                    mass_detrainment_updraft=locals.mass_detrainment_updraft,
                    mass_entrainment_u_updraft=locals.mass_entrainment_u_updraft,
                    mass_detrainment_u_updraft=locals.mass_detrainment_u_updraft,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    u_c=locals.u_c,
                    v_c=locals.v_c,
                    u_cloud_levels=locals.u_cloud_levels,
                    v_cloud_levels=locals.v_cloud_levels,
                    partition_liquid_ice=locals.partition_liquid_ice,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    vapor_excess=locals.vapor_excess,
                    t_excess=locals.t_excess,
                    add_buoyancy=locals.add_buoyancy,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # Get buoyancy of updrafts
                # NOTE test GF2020_CumulusParameterization_GetBuoyancy_2_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._get_buoyancy(
                    lcl_level=state.output.lcl_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    cloud_moist_static_energy=locals.cloud_moist_static_energy,
                    environment_moist_static_energy=locals.environment_moist_static_energy_cloud_levels,
                    environment_saturation_moist_static_energy=locals.environment_saturation_moist_static_energy_cloud_levels,
                    buoyancy=locals.buoyancy,
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # NOTE test GF2020_CumulusParameterization_GetBuoyancy_3_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._get_buoyancy(
                    lcl_level=state.output.lcl_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    cloud_moist_static_energy=locals.cloud_moist_static_energy_forced,
                    environment_moist_static_energy=locals.environment_moist_static_energy_cloud_levels_forced,
                    environment_saturation_moist_static_energy=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    buoyancy=locals.buoyancy,
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                if self.cumulus_parameterization_config.FIRST_GUESS_W:
                    # calculate in-cloud/updraft air temperature for vertical velocity
                    # NOTE test GF2020_CumulusParameterization_UpdraftTemperature_{plume}:
                    # NOTE      deep ✅
                    # NOTE      mid ✅
                    # NOTE      shallow ✅
                    self._updraft_temperature(
                        error_code=state.output.error_code,
                        updraft_column_temperature_forced=locals.updraft_column_temperature_forced,
                        cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                        geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                        cloud_total_water_after_entrainment_forced=locals.cloud_total_water_after_entrainment_forced,
                        t_cloud_levels_forced=locals.t_cloud_levels_forced,
                        plume=self.plume_dependent_constants.PLUME_INDEX,
                    )

                    # vertical velocity
                    # NOTE test GF2020_CumulusParameterization_UpdraftVerticalVelocity_{plume}:
                    # NOTE      deep ✅
                    # NOTE      mid ✅
                    # NOTE      shallow ✅
                    self._updraft_vertical_velosity(
                        vertical_velocity_3d=locals.vertical_velocity_3d,
                        vertical_velocity_2d=locals.vertical_velocity_2d,
                        convective_scale_velocity=state.input_output.convective_scale_velocity,
                        entrainment_rate=state.output.entrainment_rate,
                        detrainment_function_updraft=locals.detrainment_function_updraft,
                        geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                        t_cloud_levels_forced=locals.t_cloud_levels_forced,
                        updraft_column_temperature_forced=locals.updraft_column_temperature_forced,
                        cloud_total_water_after_entrainment_forced=locals.cloud_total_water_after_entrainment_forced,
                        cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                        vapor_forced=locals.vapor_forced,
                        updraft_lfc_level=state.output.updraft_lfc_level,
                        cloud_top_level=state.output.cloud_top_level,
                        error_code=state.output.error_code,
                        plume=self.plume_dependent_constants.PLUME_INDEX,
                    )

                # downdraft origin level
                # NOTE test GF2020_CumulusParameterization_DowndraftOriginLevel_{plume}:
                # NOTE      deep ❌ RUNS BUT DOES NOT VALIDATE. current version is brittle, likely need solver mechanics
                # NOTE      mid ❌ RUNS BUT DOES NOT VALIDATE. current version is brittle, likely need solver mechanics
                # NOTE      shallow ❌ RUNS BUT DOES NOT VALIDATE. current version is brittle, likely need solver mechanics
                self._downdraft_origin_level(
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    updraft_origin_level=state.output.updraft_origin_level,
                    downdraft_origin_level=locals.downdraft_origin_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    detrainment_start_level=locals.detrainment_start_level,
                    melting_layer=locals.melting_layer,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # downdraft normalized mass flux
                self._downdraft_normalized_mass_flux()

                # lateral mass fluxes associated with downdrafts
                # NOTE test GF2020_CumulusParameterization_DowndraftLateralMassFlux_{plume}:
                # NOTE      deep ❌ RUNS BUT DOES NOT VALIDATE
                # NOTE      mid ❌ RUNS BUT DOES NOT VALIDATE
                # NOTE      shallow ✅
                self._downdraft_lateral_mass_flux(
                    error_code=state.output.error_code,
                    downdraft_origin_level=locals.downdraft_origin_level,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_downdraft=locals.normalized_massflux_downdraft,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_updraft_forced,
                    normalized_massflux_downdraft_modified=locals.normalized_massflux_downdraft_modified,
                    detrainment_function_downdraft=locals.detrainment_function_downdraft,
                    entrainment_rate_downdraft=locals.entrainment_rate_downdraft,
                    mass_entrainment_downdraft=locals.mass_entrainment_downdraft,
                    mass_detrainment_downdraft=locals.mass_detrainment_downdraft,
                    mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    mass_entrainment_u_downdraft=locals.mass_entrainment_u_downdraft,
                    mass_detrainment_u_downdraft=locals.mass_detrainment_u_downdraft,
                    LAMBDA_DOWN=self.plume_dependent_constants.LAMBDA_DOWN,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # wet bulb temperature and moisture at downdraft origin level
                # NOTE this section does not run in the test case, and has not been implemented.
                # NOTE an error will stop execution during initalization if this would be called
                if (
                    self.cumulus_parameterization_config.USE_WETBULB
                    and self.plume_dependent_constants.PLUME_INDEX != 0
                ):
                    raise NotImplementedError(
                        "wet bulb functionality is not implemented. You should not be here,"
                        "there are multiple layers of errors that should have caught you first."
                        "If you are seeing this (at runtime), seek help."
                    )

                # downdraft moist static energy + moisture budget
                # NOTE test GF2020_CumulusParameterization_DowndraftMSEAnBuoyancy{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._downdraft_moist_static_energy_and_buoyancy(
                    error_code=state.output.error_code,
                    downdraft_origin_level=locals.downdraft_origin_level,
                    u=state.input_output.u,
                    u_cloud_levels=locals.u_cloud_levels,
                    u_c_downdraft=locals.u_c_downdraft,
                    v=state.input_output.v,
                    v_cloud_levels=locals.v_cloud_levels,
                    v_c_downdraft=locals.v_c_downdraft,
                    environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy=locals.cloud_moist_static_energy,
                    cloud_moist_static_energy_downdraft_forced=locals.cloud_moist_static_energy_downdraft_forced,
                    buoyancy_downdraft_forced=locals.d_buoyancy_downdraft_forced,
                    t_wetbulb=locals.t_wetbulb,
                    vapor_wetbulb=locals.vapor_wetbulb,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    mass_entrainment_u_downdraft=locals.mass_entrainment_u_downdraft,
                    mass_detrainment_u_downdraft=locals.mass_detrainment_u_downdraft,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # calculate moisture properties of downdraft
                # NOTE test GF2020_CumulusParameterization_DowndraftMoisture_{plume}:
                # NOTE      deep ❌ RUNS BUT DOES NOT VALIDATE
                # NOTE      mid ❌ RUNS BUT DOES NOT VALIDATE
                # NOTE      shallow ✅
                self._downdraft_moisture(
                    error_code=state.output.error_code,
                    downdraft_origin_level=locals.downdraft_origin_level,
                    t_cloud_levels_forced=locals.t_cloud_levels_forced,
                    t_wetbulb=locals.t_wetbulb,
                    vapor_forced=locals.vapor_forced,
                    vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
                    environment_saturation_mixing_ratio_cloud_levels_forced=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    cloud_total_water_after_entrainment_forced=locals.cloud_total_water_after_entrainment_forced,
                    cloud_total_water_after_entrainment_downdraft_forced=locals.cloud_total_water_after_entrainment_downdraft_forced,
                    downdraft_saturation_vapor_forced=locals.downdraft_saturation_vapor_forced,
                    vapor_wetbulb=locals.vapor_wetbulb,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    environment_moist_static_energy_forced=locals.environment_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy_downdraft_forced=locals.cloud_moist_static_energy_downdraft_forced,
                    evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    gamma_cloud_levels_forced=locals.gamma_cloud_levels_forced,
                    total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                    total_normalized_integrated_evaporate_forced=state.output.total_normalized_integrated_evaporate_forced,
                    buoyancy=locals.buoyancy,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # calculate workfunctions for updrafts
                # NOTE test GF2020_CumulusParameterization_UpdraftInitialWorkfunctions_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._updraft_initial_workfunctions(
                    error_code=state.output.error_code,
                    updraft_origin_level=state.output.updraft_origin_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_updraft=locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    d_buoyancy=locals.d_buoyancy,
                    d_buoyancy_forced=locals.d_buoyancy_forced,
                    gamma_cloud_levels=locals.gamma_cloud_levels,
                    gamma_cloud_levels_forced=locals.gamma_cloud_levels_forced,
                    t_cloud_levels=locals.t_cloud_levels,
                    t_cloud_levels_forced=locals.t_cloud_levels_forced,
                    cloud_workfunction_0=locals.cloud_workfunction_0,
                    cloud_workfunction_1=locals.cloud_workfunction_1,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # calculate CIN for updrafts
                # NOTE test GF2020_CumulusParameterization_UpdraftCIN_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._updraft_cin(
                    error_code=state.output.error_code,
                    updraft_origin_level=state.output.updraft_origin_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_updraft=locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    d_buoyancy=locals.d_buoyancy,
                    d_buoyancy_forced=locals.d_buoyancy_forced,
                    gamma_cloud_levels=locals.gamma_cloud_levels,
                    gamma_cloud_levels_forced=locals.gamma_cloud_levels_forced,
                    t_cloud_levels=locals.t_cloud_levels,
                    t_cloud_levels_forced=locals.t_cloud_levels_forced,
                    cin_0=locals.cin_0,
                    cin_1=locals.cin_1,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # trigger function: KE+CIN < 0 --> no convection
                # NOTE test GF2020_CumulusParameterization_ConvectionTrigger_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._convection_trigger(
                    error_code=state.output.error_code,
                    convective_scale_velosity=state.input_output.convective_scale_velocity,
                    cin_0=locals.cin_0,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # calculate downdraft air temperature for vertical velocitys
                # NOTE test GF2020_CumulusParameterization_DowndraftTemperature_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._downdraft_temperature(
                    state=state,
                    locals=locals,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # diurnal cycle section
                # NOTE test GF2020_CumulusParameterization_DiurnalCycle_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._diurnal_cycle(
                    error_code=state.output.error_code,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    pbl_level=state.input_output.pbl_level,
                    grid_length=state.input_output.grid_length,
                    ocean_fraction=state.input.ocean_fraction,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    t_old=state.input_output.t_old,
                    t_new=locals.t_new,
                    t_cloud_levels_forced=locals.t_cloud_levels_forced,
                    vapor_old=state.input_output.vapor_old,
                    vapor_forced=locals.vapor_forced,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    vertical_velocity_2d=locals.vertical_velocity_2d,
                    cape_removal_time_scale=locals.cape_removal_time_scale,
                    cape_removal_time_scale_from_state=state.output.cape_removal_time_scale,
                    pbl_time_scale=locals.pbl_time_scale,
                    pbl_time_scale_from_state=state.output.pbl_time_scale,
                    cloud_work_function_1_pbl=locals.cloud_work_function_1_pbl,
                    cloud_work_function_1_fa=locals.cloud_work_function_1_fa,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # Trigger function based on Xie et al 2019
                # NOTE not implemented, does not run with test config
                self._Xie_trigger_function()

                # determine downdraft strength in terms of windshear
                # NOTE test GF2020_CumulusParameterization_DowndraftWindshear_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ❌ one field, one point (0.17%), 4 ULP
                # NOTE      shallow ✅
                self._downdraft_windshear(
                    error_code=state.output.error_code,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    geopotential_height_forced=state.input_output.geopotential_height_forced,
                    p_forced=state.input_output.p_forced,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    ccn=state.input_output.ccn,
                    psum=locals.psum,
                    psumh=locals.psumh,
                    total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                    total_normalized_integrated_evaporate_forced=state.output.total_normalized_integrated_evaporate_forced,
                    scale_dependence_factor_downdraft=locals.scale_dependence_factor_downdraft,
                    epsilon=locals.epsilon,
                    epsilon_min=locals.epsilon_min,
                    epsilon_max=locals.epsilon_max,
                    epsilon_computed=locals.epsilon_computed,
                    epsilon_forced=state.output.epsilon_forced,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # get the environmental mass flux
                # NOTE test GF2020_CumulusParameterization_EnvironmentMassFlux_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._environment_mass_flux(
                    error_code=state.output.error_code,
                    epsilon_forced=state.output.epsilon_forced,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    environment_massflux=locals.environment_massflux,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # check mass conservation
                # NOTE This code runs in the Fortran and only has one output: totmas (total mass).
                # totmas has only one use: a conditional log write if total mass is above a 1e-6.
                # This conditional also has a disabled fatal error call.
                # Since the only consequential outcome is disabled, and this port has thus far not
                # implemented other log writes, this code not been implemented.
                # If totmas is needed in the future, or this fatal call is reimplemented,
                # this code will be revisited
                self._mass_conservation()

                # change per unit mass that a model cloud would modify the environment
                # NOTE test GF2020_CumulusParameterization_VerticalDiscretization_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._vertical_discretization(
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    environment_massflux=locals.environment_massflux,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    c1d=locals.c1d,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    u_cloud_levels=locals.u_cloud_levels,
                    v_cloud_levels=locals.v_cloud_levels,
                    u_c=locals.u_c,
                    v_c=locals.v_c,
                    u_c_downdraft=locals.u_c_downdraft,
                    v_c_downdraft=locals.v_c_downdraft,
                    cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                    cloud_moist_static_energy_downdraft_forced=locals.cloud_moist_static_energy_downdraft_forced,
                    environment_moist_static_energy_cloud_levels_forced=locals.environment_moist_static_energy_cloud_levels_forced,
                    vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
                    cloud_total_water_after_entrainment_forced=locals.cloud_total_water_after_entrainment_forced,
                    cloud_total_water_after_entrainment_downdraft_forced=locals.cloud_total_water_after_entrainment_downdraft_forced,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                    melting=locals.melting,
                    partition_liquid_ice=locals.partition_liquid_ice,
                    epsilon_forced=state.output.epsilon_forced,
                    d_buoyancy_downdraft_forced=locals.d_buoyancy_downdraft_forced,
                    del_u_cloud_ensemble=locals.del_u_cloud_ensemble,
                    del_v_cloud_ensemble=locals.del_v_cloud_ensemble,
                    del_moist_static_energy_cloud_ensemble=locals.del_moist_static_energy_cloud_ensemble,
                    del_t_cloud_ensemble=locals.del_t_cloud_ensemble,
                    del_vapor_cloud_ensemble=locals.del_vapor_cloud_ensemble,
                    del_cloud_liquid_cloud_ensemble=locals.del_cloud_liquid_cloud_ensemble,
                    del_buoyancy_cloud_ensemble=locals.del_buoyancy_cloud_ensemble,
                    t_tendency_from_environmental_subsidence=locals.t_tendency_from_environmental_subsidence,
                    moist_static_energy_tendency_from_environmental_subsidence=locals.moist_static_energy_tendency_from_environmental_subsidence,
                    vapor_tendency_from_environmental_subsidence=locals.vapor_tendency_from_environmental_subsidence,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

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
