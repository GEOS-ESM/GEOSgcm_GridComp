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
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MAXENS1, MAXENS2, MAXENS3
from pyMoist.convection.GF_2020.cumulus_parameterization.environment import (
    environment_conditions,
    environment_cloud_levels,
    environment_mass_flux,
    modify_environment_profiles,
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
    get_precip_fluxes,
    PrecipFactor,
    rain_evaporation_below_cloud_base,
    cloud_dissapation,
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
from pyMoist.convection.GF_2020.cumulus_parameterization.convective_tracers import (
    ColdPoolParameterization,
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
    StaticControl,
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
    UpdraftWorkfunctions,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.triggers import (
    convection_trigger,
    XieTriggerFunction,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.downdraft import (
    DowndraftOriginLevel,
    downdraft_mass_flux,
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
from pyMoist.convection.GF_2020.cumulus_parameterization.smoothing import smooth_tendencies
from pyMoist.convection.GF_2020.cumulus_parameterization.large_scale_forcing import (
    LargeScaleForcing,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.kinetic_energy_to_heating import (
    kinetic_energy_to_heating,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.feedback.feedback import Feedback
from pyMoist.convection.GF_2020.cumulus_parameterization.prepare_output import (
    total_evaporation_flux,
    LightningFlashDensity,
    deep_precipitation_output,
    tracer_output,
    prepare_output,
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

        # initalize the local fields needed within this class
        self.locals = GF2020CumulusParameterizationLocals.zeros(
            quantity_factory,
            data_dimensions={
                "ensemble_1": MAXENS1,
                "ensemble_2": MAXENS2,
                "ensemble_3": MAXENS3,
                "ensemble_members": MAXENS1 * MAXENS2 * MAXENS3,
                "tracers": config.NUMBER_OF_TRACERS,
            },
        )

        # initalize all the subclasses
        self._setup = Setup(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._environment_conditions = stencil_factory.from_dims_halo(
            func=environment_conditions,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "SATURATION_CALCULATION_CHOICE": cumulus_parameterization_config.SATURATION_CALCULATION_CHOICE
            },
        )

        self._sounding = Sounding(cumulus_parameterization_config=cumulus_parameterization_config)

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

        self._cold_pool_parameterization = ColdPoolParameterization(
            cumulus_parameterization_config=cumulus_parameterization_config
        )

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
            externals={
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                "MIN_ENTRAINMENT_RATE": cumulus_parameterization_config.MIN_ENTRAINMENT_RATE,
            },
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

        self._updraft_vertical_velocity = stencil_factory.from_dims_halo(
            func=updraft_vertical_velocity,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF},
        )

        self._downdraft_origin_level = DowndraftOriginLevel(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._downdraft_mass_flux = stencil_factory.from_dims_halo(
            func=downdraft_mass_flux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF},
        )

        self._downdraft_lateral_mass_flux = stencil_factory.from_dims_halo(
            func=downdraft_lateral_massflux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._downdraft_wet_bulb = DowndraftWetBlub(
            cumulus_parameterization_config=cumulus_parameterization_config
        )

        self._downdraft_moist_static_energy_and_buoyancy = stencil_factory.from_dims_halo(
            func=downdraft_moist_static_energy_and_buoyancy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_WETBULB": cumulus_parameterization_config.USE_WETBULB,
                "PRESSURE_GRADIENT_CONSTANT": cumulus_parameterization_config.PRESSURE_GRADIENT_CONSTANT,
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

        self._Xie_trigger_function = XieTriggerFunction(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

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

        self._environmental_subsidence = EnvironmentalSubsidence(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._smooth_tendencies = stencil_factory.from_dims_halo(
            func=smooth_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"USE_SMOOTH_TENDENCIES": cumulus_parameterization_config.USE_SMOOTH_TENDENCIES},
        )

        self._modify_environment_profiles = stencil_factory.from_dims_halo(
            func=modify_environment_profiles,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "COUPLE_MICROPHYSICS": cumulus_parameterization_config.COUPLE_MICROPHYSICS,
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
            },
        )
        self._static_control = StaticControl(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._updraft_workfunctions = UpdraftWorkfunctions(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._large_scale_forcing = LargeScaleForcing(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._kinetic_energy_to_heating = stencil_factory.from_dims_halo(
            func=kinetic_energy_to_heating,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._feedback = Feedback()

        self._precipitation_flux = stencil_factory.from_dims_halo(
            func=get_precip_fluxes,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._rain_evaporation_below_cloud_base = stencil_factory.from_dims_halo(
            func=rain_evaporation_below_cloud_base,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._cloud_dissapation = stencil_factory.from_dims_halo(
            func=cloud_dissapation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "USE_CLOUD_DISSIPATION": cumulus_parameterization_config.USE_CLOUD_DISSIPATION,
                "DTIME": cumulus_parameterization_config.DTIME,
                "COUPLE_MICROPHYSICS": cumulus_parameterization_config.COUPLE_MICROPHYSICS,
            },
        )

        self._total_evaporation_flux = stencil_factory.from_dims_halo(
            func=total_evaporation_flux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._lightning_flash_density = LightningFlashDensity(
            cumulus_parameterization_config=cumulus_parameterization_config
        )

        self._deep_precipitation_output = stencil_factory.from_dims_halo(
            func=deep_precipitation_output,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._tracer_output = stencil_factory.from_dims_halo(
            func=tracer_output,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._prepare_output = stencil_factory.from_dims_halo(
            func=prepare_output,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._atmospheric_composition = AtmosphericComposition(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
        )

        self._gate_sounding = GATESounding(cumulus_parameterization_config=cumulus_parameterization_config)

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        saturation_tables: SaturationVaporPressureTable,
    ):
        if self.cumulus_parameterization_config.PLUME_ORDER == 0:
            plume_types = ["shallow", "mid", "deep"]
        elif self.cumulus_parameterization_config.PLUME_ORDER == 1:
            plume_types = ["shallow", "deep", "mid"]
        else:
            raise NotImplementedError("plume order not impelemented")

        for plume in plume_types:
            # setup constants for the current plume, reset necessary fields, prefill necessary fields
            # NOTE test GF2020_CumulusParameterization_Setup_{plume}:
            # NOTE      deep ✅
            # NOTE      mid ✅
            # NOTE      shallow ❌ (cannot serialize output data because the plume is disabled)
            self._setup(
                error_code=state.output.error_code,
                error_code_2=self.locals.error_code_2,
                error_code_3=self.locals.error_code_3,
                maximum_updraft_origin_level=self.locals.maximum_updraft_origin_level,
                kstabm=self.locals.kstabm,
                t_excess=state.input.t_excess,
                t_excess_local=self.locals.t_excess,
                vapor_excess=state.input.vapor_excess,
                vapor_excess_local=self.locals.vapor_excess,
                ocean_fraction=state.input.ocean_fraction,
                ocean_fraction_local=self.locals.ocean_fraction,
                t_old=state.input_output.t_old,
                t_new=self.locals.t_new,
                t_new_pbl=self.locals.t_new_pbl,
                vapor_old=state.input_output.vapor_old,
                vapor_forced=self.locals.vapor_forced,
                vapor_forced_pbl=self.locals.vapor_forced_pbl,
                downdraft_saturation_vapor_forced=self.locals.downdraft_saturation_vapor_forced,
                grid_scale_forcing_t=state.input.grid_scale_forcing_t,
                grid_scale_forcing_vapor=state.input.grid_scale_forcing_vapor,
                subgrid_scale_forcing_t=state.input.subgrid_scale_forcing_t,
                subgrid_scale_forcing_vapor=state.input.subgrid_scale_forcing_vapor,
                dmoist_static_energydt=self.locals.dmoist_static_energydt,
                cloud_moist_static_energy_downdraft_forced=self.locals.cloud_moist_static_energy_downdraft_forced,
                cloud_moist_static_energy_forced_transported=self.locals.cloud_moist_static_energy_forced_transported,
                cap_max=self.locals.cap_max,
                cap_max_increment=self.locals.cap_max_increment,
                geopotential_height=state.input_output.geopotential_height_forced,
                geopotential_height_local=self.locals.geopotential_height,
                geopotential_height_modified_local=self.locals.geopotential_height_modified,
                cloud_workfunction_0=self.locals.cloud_workfunction_0,
                cloud_workfunction_1=self.locals.cloud_workfunction_1,
                cloud_workfunction_2=self.locals.cloud_workfunction_2,
                cloud_workfunction_3=self.locals.cloud_workfunction_3,
                cloud_workfunction_0_pbl=self.locals.cloud_workfunction_0_pbl,
                cloud_workfunction_1_pbl=self.locals.cloud_workfunction_1_pbl,
                cloud_workfunction_1_fa=self.locals.cloud_workfunction_1_fa,
                cin_1=self.locals.cin_1,
                k_x_modified=self.locals.k_x_modified,
                epsilon_forced=state.output.epsilon_forced,
                epsilon_local=self.locals.epsilon,
                epsilon_min=self.locals.epsilon_min,
                epsilon_max=self.locals.epsilon_max,
                pbl_time_scale=self.locals.pbl_time_scale,
                t_wetbulb=self.locals.t_wetbulb,
                vapor_wetbulb=self.locals.vapor_wetbulb,
                cape_removal_time_scale=self.locals.cape_removal_time_scale,
                f_dicycle_modified=self.locals.f_dicycle_modified,
                add_buoyancy=self.locals.add_buoyancy,
                scale_dependence_factor=state.output.scale_dependence_factor,
                scale_dependence_factor_downdraft=self.locals.scale_dependence_factor_downdraft,
                c1d=self.locals.c1d,
                evaporation_below_cloud_base=self.locals.evaporation_below_cloud_base,
                mass_flux_ensemble=self.locals.mass_flux_ensemble,
                precipitation_ensemble=self.locals.precipitation_ensemble,
                precip=state.output.precip,
                lightning_density=state.output.lightning_density,
                seed_convection=state.input.seed_convection,
                grid_length=state.input_output.grid_length,
                random_number=self.locals.random_number,
                lateral_entrainment_rate=state.input.lateral_entrainment_rate,
                entrainment_rate=state.output.entrainment_rate,
                detrainment_function_updraft=self.locals.detrainment_function_updraft,
                arbitrary_numerical_parameter=self.locals.arbitrary_numerical_parameter,
                plume_dependent_constants=self.plume_dependent_constants,
                plume=plume,
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
                    moist_static_energy=self.locals.environment_moist_static_energy,
                    saturation_moist_static_energy=self.locals.environment_saturation_moist_static_energy,
                    saturation_mixing_ratio=self.locals.environment_saturation_mixing_ratio,
                    geopotential_height=self.locals.geopotential_height,
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
                    t=self.locals.t_new,
                    vapor=self.locals.vapor_forced,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    moist_static_energy=self.locals.environment_moist_static_energy_forced,
                    saturation_moist_static_energy=self.locals.environment_saturation_moist_static_energy_forced,
                    saturation_mixing_ratio=self.locals.environment_saturation_mixing_ratio_forced,
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
                    p_cloud_levels=self.locals.p_cloud_levels,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    geopotential_height=self.locals.geopotential_height,
                    geopotential_height_cloud_levels=self.locals.geopotential_height_cloud_levels,
                    t=state.input_output.t_old,
                    t_surface=state.input_output.t_surface,
                    t_cloud_levels=self.locals.t_cloud_levels,
                    vapor=state.input_output.vapor_old,
                    vapor_cloud_levels=self.locals.vapor_cloud_levels,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    u_cloud_levels=self.locals.u_cloud_levels,
                    v_cloud_levels=self.locals.v_cloud_levels,
                    environment_saturation_mixing_ratio=self.locals.environment_saturation_mixing_ratio,
                    environment_saturation_mixing_ratio_cloud_levels=self.locals.environment_saturation_mixing_ratio_cloud_levels,
                    environment_moist_static_energy=self.locals.environment_moist_static_energy,
                    environment_moist_static_energy_cloud_levels=self.locals.environment_moist_static_energy_cloud_levels,
                    environment_saturation_moist_static_energy=self.locals.environment_saturation_moist_static_energy,
                    environment_saturation_moist_static_energy_cloud_levels=self.locals.environment_saturation_moist_static_energy_cloud_levels,
                    gamma_cloud_levels=self.locals.gamma_cloud_levels,
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
                    geopotential_height_cloud_levels=self.locals.geopotential_height_cloud_levels_forced,
                    t=self.locals.t_new,
                    t_surface=state.input_output.t_surface,
                    t_cloud_levels=self.locals.t_cloud_levels_forced,
                    vapor=self.locals.vapor_forced,
                    vapor_cloud_levels=self.locals.vapor_cloud_levels_forced,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    u_cloud_levels=self.locals.u_cloud_levels,
                    v_cloud_levels=self.locals.v_cloud_levels,
                    environment_saturation_mixing_ratio=self.locals.environment_saturation_mixing_ratio_forced,
                    environment_saturation_mixing_ratio_cloud_levels=self.locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    environment_moist_static_energy=self.locals.environment_moist_static_energy_forced,
                    environment_moist_static_energy_cloud_levels=self.locals.environment_moist_static_energy_cloud_levels_forced,
                    environment_saturation_moist_static_energy=self.locals.environment_saturation_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    gamma_cloud_levels=self.locals.gamma_cloud_levels_forced,
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
                    geopotential_height=self.locals.geopotential_height_cloud_levels_forced,
                    error_code=state.output.error_code,
                    air_density=self.locals.hydrostatic_air_density,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # partition between liq/ice cloud contents
                # NOTE test GF2020_CumulusParameterization_PartitionLiquidIceAndGetLevels_{plume}:
                # NOTE      deep ❌ 1 var failing (500 ULP) -> ice_fraction function causing issues
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._partition_liquid_ice(
                    t=self.locals.t_new,
                    p=state.output.p_cloud_levels_forced,
                    geopotential_height=self.locals.geopotential_height_cloud_levels_forced,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    surface_type=state.input.surface_type,
                    convection_fraction=state.input.convection_fraction,
                    error_code=state.output.error_code,
                    melting_layer=self.locals.melting_layer,
                    part_liquid_ice=self.locals.partition_liquid_ice,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                self._find_maximum_updraft_origin_level(
                    geopotential_height=self.locals.geopotential_height_cloud_levels_forced,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    error_code=state.output.error_code,
                    maximum_updraft_origin_level=self.locals.maximum_updraft_origin_level,
                    MAX_UPDRAFT_ORIGIN_HEIGHT=self.plume_dependent_constants.MAX_UPDRAFT_ORIGIN_HEIGHT,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                self._find_detrainmet_start_level(
                    geopotential_height=self.locals.geopotential_height_cloud_levels_forced,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    error_code=state.output.error_code,
                    detrainment_start_level=self.locals.detrainment_start_level,
                    DETRAINMENT_CRITICAL_DEPTH=self.plume_dependent_constants.DETRAINMENT_CRITICAL_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # determine level with highest moist static energy content (k_max_mse)
                # NOTE test GF2020_CumulusParameterization_HighestMoistStaticEnergyLevel_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._find_highest_moist_static_energy_level(
                    moist_static_energy=self.locals.environment_moist_static_energy_cloud_levels_forced,
                    error_code=state.output.error_code,
                    maximum_updraft_origin_level=self.locals.maximum_updraft_origin_level,
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
                    p_cloud_levels=self.locals.p_cloud_levels,
                    t_excess=self.locals.t_excess,
                    t_cloud_levels_forced=self.locals.t_cloud_levels,
                    t_perturbation=state.output.t_perturbation,
                    vapor_excess=self.locals.vapor_excess,
                    vapor_cloud_levels_forced=self.locals.vapor_cloud_levels,
                    omega=state.input_output.omega,
                    air_density=state.input_output.air_density,
                    geopotential_height_cloud_levels=self.locals.geopotential_height_cloud_levels,
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
                    t_excess=self.locals.t_excess,
                    vapor_excess=self.locals.vapor_excess,
                    add_buoyancy=self.locals.add_buoyancy,
                    ocean_fraction=state.input.ocean_fraction,
                    updraft_origin_level=state.output.updraft_origin_level,
                    p=state.input_output.p_forced,
                    environmenet_moist_static_energy=self.locals.environment_moist_static_energy_cloud_levels,
                    environmenet_moist_static_energy_forced=self.locals.environment_moist_static_energy_cloud_levels_forced,
                    t_perturbation=state.output.t_perturbation,
                    moist_static_energy_origin_level=self.locals.moist_static_energy_origin_level,
                    moist_static_energy_origin_level_forced=self.locals.moist_static_energy_origin_level_forced,
                    AVERAGE_LAYER_DEPTH=self.plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # determine the vertical entrainment/detrainment rates
                # NOTE test GF2020_CumulusParameterization_EntrainmentRates_{plume}:
                # NOTE      deep ❌ (worst fail 2.41% - max 4 ULP)
                # NOTE      mid ❌ (worst fail 0.88% - max 2 ULP)
                # NOTE      shallow ✅
                self._entrainment_rates(
                    vapor=self.locals.vapor_cloud_levels_forced,
                    environment_saturation_mixing_ratio=self.locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    lcl_level=state.output.lcl_level,
                    error_code=state.output.error_code,
                    entrainment_rate=state.output.entrainment_rate,
                    detrainment_function_updraft=self.locals.detrainment_function_updraft,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # determine level of convective cloud base
                # NOTE test GF2020_CumulusParameterization_ConvectiveCloudBaseLevel_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._set_start_level(
                    lcl_level=state.output.lcl_level,
                    start_level=self.locals.start_level,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                self._get_convective_cloud_base_level(
                    error_code=state.output.error_code,
                    lcl_level=state.output.lcl_level,
                    cloud_moist_static_energy_forced_transported=self.locals.cloud_moist_static_energy_forced_transported,
                    cap_max=self.locals.cap_max,
                    updraft_origin_level=state.output.updraft_origin_level,
                    start_level=self.locals.start_level,
                    moist_static_energy_origin_level_forced=self.locals.moist_static_energy_origin_level_forced,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    maximum_updraft_origin_level=self.locals.maximum_updraft_origin_level,
                    negative_buoyancy_depth=self.locals.negative_buoyancy_depth,
                    frh_lfc=self.locals.frh_lfc,
                    geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                    entrainment_rate=state.output.entrainment_rate,
                    environment_moist_static_energy_forced=self.locals.environment_moist_static_energy_forced,
                    environment_moist_static_energy_cloud_levels_forced=self.locals.environment_moist_static_energy_cloud_levels_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    t_excess=self.locals.t_excess,
                    vapor_excess=self.locals.vapor_excess,
                    add_buoyancy=self.locals.add_buoyancy,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    vapor_forced=self.locals.vapor_forced,
                    environment_saturation_mixing_ratio_forced=self.locals.environment_saturation_mixing_ratio_forced,
                    ocean_fraction=self.locals.ocean_fraction,
                    cap_max_increment=self.locals.cap_max_increment,
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
                    entrainment_rate_downdraft=self.locals.entrainment_rate_downdraft,
                    detrainment_function_downdraft=self.locals.detrainment_function_downdraft,
                    scale_dependence_factor_downdraft=self.locals.scale_dependence_factor_downdraft,
                    plume_entrainment_rate=self.plume_dependent_constants.ENTRAINMENT_RATE,
                )

                # update unforced & forced moist static energy
                # NOTE test GF2020_CumulusParameterization_ParcelMoistStaticEnergy_2_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._parcel_moist_static_energy(
                    error_code=state.output.error_code,
                    t_excess=self.locals.t_excess,
                    vapor_excess=self.locals.vapor_excess,
                    add_buoyancy=self.locals.add_buoyancy,
                    ocean_fraction=state.input.ocean_fraction,
                    updraft_origin_level=state.output.updraft_origin_level,
                    p=state.input_output.p_forced,
                    environmenet_moist_static_energy=self.locals.environment_moist_static_energy_cloud_levels,
                    environmenet_moist_static_energy_forced=self.locals.environment_moist_static_energy_cloud_levels_forced,
                    t_perturbation=state.output.t_perturbation,
                    moist_static_energy_origin_level=self.locals.moist_static_energy_origin_level,
                    moist_static_energy_origin_level_forced=self.locals.moist_static_energy_origin_level_forced,
                    AVERAGE_LAYER_DEPTH=self.plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # increase detrainment in stable layers
                # NOTE test GF2020_CumulusParameterization_StableDetrainment_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._unknown_find_level(
                    array=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    start_index=state.output.updraft_lfc_level,
                    end_index=self.locals.kstabm,
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
                    moist_static_energy=self.locals.environment_moist_static_energy_forced,
                    saturation_moist_static_energy=self.locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    moist_static_energy_origin_level=self.locals.moist_static_energy_origin_level_forced,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    geopotential_height=self.locals.geopotential_height_cloud_levels_forced,
                    cloud_moist_static_energy=self.locals.cloud_moist_static_energy_forced_transported,
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                self._cloud_top_checks(
                    cloud_top_level=state.output.cloud_top_level,
                    p=state.output.p_cloud_levels_forced,
                    geopotential_height=self.locals.geopotential_height_cloud_levels,
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
                    error_code=state.output.error_code,
                    updraft_origin_level=state.output.updraft_origin_level,
                    cloud_top_level=state.output.cloud_top_level,
                    pbl_level=state.input_output.pbl_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    lcl_level=state.output.lcl_level,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    p_surface=state.input_output.p_surface,
                    ocean_fraction=state.input.ocean_fraction,
                    normalized_massflux_updraft=self.locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    normalized_massflux_updraft_modified=self.locals.normalized_massflux_updraft_modified,
                    random_number=self.locals.random_number,
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
                    geopotential_height=self.locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_updraft=state.output.normalized_massflux_updraft_forced,
                    detrainment_function_updraft=self.locals.detrainment_function_updraft,
                    entrainment_rate=state.output.entrainment_rate,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_entrainment_updraft=self.locals.mass_entrainment_updraft,
                    mass_detrainment_updraft=self.locals.mass_detrainment_updraft,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    updraft_origin_level=state.output.updraft_origin_level,
                    pbl_level=state.input_output.pbl_level,
                    mass_entrainment_u_updraft=self.locals.mass_entrainment_u_updraft,
                    mass_detrainment_u_updraft=self.locals.mass_detrainment_u_updraft,
                    LAMBDA_DEEP=self.plume_dependent_constants.LAMBDA_DEEP,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                self._compute_uc_vc(
                    u_c=self.locals.u_c,
                    v_c=self.locals.v_c,
                    cloud_moist_static_energy=self.locals.cloud_moist_static_energy,
                    cloud_moist_static_energy_forced=self.locals.cloud_moist_static_energy_forced,
                    error_code=state.output.error_code,
                    start_level=self.locals.start_level,
                    moist_static_energy_origin_level=self.locals.moist_static_energy_origin_level,
                    moist_static_energy_origin_level_forced=self.locals.moist_static_energy_origin_level_forced,
                    u_cloud_levels=self.locals.u_cloud_levels,
                    v_cloud_levels=self.locals.v_cloud_levels,
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
                    start_level=self.locals.start_level,
                    cloud_top_level=state.output.cloud_top_level,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                    normalized_massflux_updraft=self.locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    environment_moist_static_energy_forced=self.locals.environment_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy_forced=self.locals.cloud_moist_static_energy_forced,
                    vapor_excess=self.locals.vapor_excess,
                    t_excess=self.locals.t_excess,
                    add_buoyancy=self.locals.add_buoyancy,
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
                    cloud_moist_static_energy=self.locals.cloud_moist_static_energy_forced,
                    environment_moist_static_energy=self.locals.environment_moist_static_energy_cloud_levels_forced,
                    environment_saturation_moist_static_energy=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    buoyancy=self.locals.buoyancy,
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
                        locals=self.locals,
                        plume_dependent_constants=self.plume_dependent_constants,
                    )

                # calculate moisture properties of updraft
                # NOTE test GF2020_CumulusParameterization_UpdraftMoisture_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._updraft_moisture_profile(
                    start_level=self.locals.start_level,
                    error_code=state.output.error_code,
                    geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                    cloud_total_water_after_entrainment_forced=self.locals.cloud_total_water_after_entrainment_forced,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                    cloud_moist_static_energy_forced=self.locals.cloud_moist_static_energy_forced,
                    updraft_column_temperature_forced=self.locals.updraft_column_temperature_forced,
                    ocean_fraction=state.input.ocean_fraction,
                    convection_fraction=state.input.convection_fraction,
                    surface_type=state.input.surface_type,
                    p_forced=state.input_output.p_forced,
                    cloud_top_level=state.output.cloud_top_level,
                    d_buoyancy_forced=self.locals.d_buoyancy_forced,
                    cloud_liquid_before_rain_forced=self.locals.cloud_liquid_before_rain_forced,
                    t_cloud_levels=self.locals.t_cloud_levels,
                    vapor_forced=self.locals.vapor_forced,
                    gamma_cloud_levels_forced=self.locals.gamma_cloud_levels_forced,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    environment_saturation_mixing_ratio_cloud_levels_forced=self.locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    updraft_origin_level=state.output.updraft_origin_level,
                    vapor_cloud_levels_forced=self.locals.vapor_cloud_levels_forced,
                    vapor_excess=self.locals.vapor_excess,
                    ccn=state.input_output.ccn,
                    mass_entrainment_updraft=self.locals.mass_entrainment_updraft,
                    mass_detrainment_updraft=self.locals.mass_detrainment_updraft,
                    psum=self.locals.psum,
                    psumh=self.locals.psumh,
                    c1d=self.locals.c1d,
                    add_buoyancy=self.locals.add_buoyancy,
                    vertical_velocity_3d=self.locals.vertical_velocity_3d,
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
                    melting_layer=self.locals.melting_layer,
                    partition_liquid_ice=self.locals.partition_liquid_ice,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    melting=self.locals.melting,
                )

                # updraft moist static energy + momentum budget
                # NOTE test GF2020_CumulusParameterization_UpdraftMoistStaticEnergyAndMomentumBudget_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._updraft_moist_static_energy_and_momentum_budget(
                    error_code=state.output.error_code,
                    start_level=self.locals.start_level,
                    cloud_top_level=state.output.cloud_top_level,
                    p_forced=state.input_output.p_forced,
                    environment_moist_static_energy=self.locals.environment_moist_static_energy,
                    environment_moist_static_energy_forced=self.locals.environment_moist_static_energy_forced,
                    environment_moist_static_energy_cloud_levels=self.locals.environment_moist_static_energy_cloud_levels,
                    environment_moist_static_energy_cloud_levels_forced=self.locals.environment_moist_static_energy_cloud_levels_forced,
                    environment_saturation_moist_static_energy_cloud_levels=self.locals.environment_saturation_moist_static_energy_cloud_levels,
                    environment_saturation_moist_static_energy_cloud_levels_forced=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy=self.locals.cloud_moist_static_energy,
                    cloud_moist_static_energy_forced=self.locals.cloud_moist_static_energy_forced,
                    normalized_massflux_updraft=self.locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    mass_entrainment_updraft=self.locals.mass_entrainment_updraft,
                    mass_detrainment_updraft=self.locals.mass_detrainment_updraft,
                    mass_entrainment_u_updraft=self.locals.mass_entrainment_u_updraft,
                    mass_detrainment_u_updraft=self.locals.mass_detrainment_u_updraft,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    u_c=self.locals.u_c,
                    v_c=self.locals.v_c,
                    u_cloud_levels=self.locals.u_cloud_levels,
                    v_cloud_levels=self.locals.v_cloud_levels,
                    partition_liquid_ice=self.locals.partition_liquid_ice,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    vapor_excess=self.locals.vapor_excess,
                    t_excess=self.locals.t_excess,
                    add_buoyancy=self.locals.add_buoyancy,
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
                    cloud_moist_static_energy=self.locals.cloud_moist_static_energy,
                    environment_moist_static_energy=self.locals.environment_moist_static_energy_cloud_levels,
                    environment_saturation_moist_static_energy=self.locals.environment_saturation_moist_static_energy_cloud_levels,
                    buoyancy=self.locals.buoyancy,
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
                    cloud_moist_static_energy=self.locals.cloud_moist_static_energy_forced,
                    environment_moist_static_energy=self.locals.environment_moist_static_energy_cloud_levels_forced,
                    environment_saturation_moist_static_energy=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    buoyancy=self.locals.buoyancy,
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
                        updraft_column_temperature_forced=self.locals.updraft_column_temperature_forced,
                        cloud_moist_static_energy_forced=self.locals.cloud_moist_static_energy_forced,
                        geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                        cloud_total_water_after_entrainment_forced=self.locals.cloud_total_water_after_entrainment_forced,
                        t_cloud_levels_forced=self.locals.t_cloud_levels_forced,
                        plume=self.plume_dependent_constants.PLUME_INDEX,
                    )

                    # vertical velocity
                    # NOTE test GF2020_CumulusParameterization_UpdraftVerticalVelocity_{plume}:
                    # NOTE      deep ✅
                    # NOTE      mid ✅
                    # NOTE      shallow ✅
                    self._updraft_vertical_velocity(
                        vertical_velocity_3d=self.locals.vertical_velocity_3d,
                        vertical_velocity_2d=self.locals.vertical_velocity_2d,
                        convective_scale_velocity=state.input_output.convective_scale_velocity,
                        entrainment_rate=state.output.entrainment_rate,
                        detrainment_function_updraft=self.locals.detrainment_function_updraft,
                        geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                        t_cloud_levels_forced=self.locals.t_cloud_levels_forced,
                        updraft_column_temperature_forced=self.locals.updraft_column_temperature_forced,
                        cloud_total_water_after_entrainment_forced=self.locals.cloud_total_water_after_entrainment_forced,
                        cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                        vapor_forced=self.locals.vapor_forced,
                        updraft_lfc_level=state.output.updraft_lfc_level,
                        cloud_top_level=state.output.cloud_top_level,
                        error_code=state.output.error_code,
                        plume=self.plume_dependent_constants.PLUME_INDEX,
                    )

                # downdraft origin level
                # NOTE test GF2020_CumulusParameterization_DowndraftOriginLevel_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._downdraft_origin_level(
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    environment_saturation_moist_static_energy_cloud_levels_forced=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    updraft_origin_level=state.output.updraft_origin_level,
                    downdraft_origin_level=state.output.downdraft_origin_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    detrainment_start_level=self.locals.detrainment_start_level,
                    melting_layer=self.locals.melting_layer,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # downdraft normalized mass flux
                # NOTE test GF2020_CumulusParameterization_DowndraftMassFlux_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._downdraft_mass_flux(
                    error_code=state.output.error_code,
                    detrainment_start_level=self.locals.detrainment_start_level,
                    downdraft_origin_level=state.output.downdraft_origin_level,
                    pbl_level=state.input_output.pbl_level,
                    updraft_origin_level=state.output.updraft_origin_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    lcl_level=state.output.lcl_level,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    p_surface=state.input_output.p_surface,
                    normalized_massflux_downdraft=self.locals.normalized_massflux_downdraft,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    ocean_fraction=state.input.ocean_fraction,
                    random_number=self.locals.random_number,
                    DOWNDRAFT_MAX_HEIGHT_LAND=self.plume_dependent_constants.DOWNDRAFT_MAX_HEIGHT_LAND,
                    DOWNDRAFT_MAX_HEIGHT_OCEAN=self.plume_dependent_constants.DOWNDRAFT_MAX_HEIGHT_OCEAN,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # lateral mass fluxes associated with downdrafts
                # NOTE test GF2020_CumulusParameterization_DowndraftLateralMassFlux_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._downdraft_lateral_mass_flux(
                    error_code=state.output.error_code,
                    downdraft_origin_level=state.output.downdraft_origin_level,
                    geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_downdraft=self.locals.normalized_massflux_downdraft,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_updraft_forced,
                    normalized_massflux_downdraft_modified=self.locals.normalized_massflux_downdraft_modified,
                    detrainment_function_downdraft=self.locals.detrainment_function_downdraft,
                    entrainment_rate_downdraft=self.locals.entrainment_rate_downdraft,
                    mass_entrainment_downdraft=self.locals.mass_entrainment_downdraft,
                    mass_detrainment_downdraft=self.locals.mass_detrainment_downdraft,
                    mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    mass_entrainment_u_downdraft=self.locals.mass_entrainment_u_downdraft,
                    mass_detrainment_u_downdraft=self.locals.mass_detrainment_u_downdraft,
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
                    downdraft_origin_level=state.output.downdraft_origin_level,
                    u=state.input_output.u,
                    u_cloud_levels=self.locals.u_cloud_levels,
                    u_c_downdraft=self.locals.u_c_downdraft,
                    v=state.input_output.v,
                    v_cloud_levels=self.locals.v_cloud_levels,
                    v_c_downdraft=self.locals.v_c_downdraft,
                    environment_moist_static_energy_forced=self.locals.environment_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy=self.locals.cloud_moist_static_energy,
                    cloud_moist_static_energy_downdraft_forced=self.locals.cloud_moist_static_energy_downdraft_forced,
                    buoyancy_downdraft_forced=self.locals.d_buoyancy_downdraft_forced,
                    t_wetbulb=self.locals.t_wetbulb,
                    vapor_wetbulb=self.locals.vapor_wetbulb,
                    geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    mass_entrainment_u_downdraft=self.locals.mass_entrainment_u_downdraft,
                    mass_detrainment_u_downdraft=self.locals.mass_detrainment_u_downdraft,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # calculate moisture properties of downdraft
                # NOTE test GF2020_CumulusParameterization_DowndraftMoisture_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._downdraft_moisture(
                    error_code=state.output.error_code,
                    downdraft_origin_level=state.output.downdraft_origin_level,
                    t_cloud_levels_forced=self.locals.t_cloud_levels_forced,
                    t_wetbulb=self.locals.t_wetbulb,
                    vapor_forced=self.locals.vapor_forced,
                    vapor_cloud_levels_forced=self.locals.vapor_cloud_levels_forced,
                    environment_saturation_mixing_ratio_cloud_levels_forced=self.locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    cloud_total_water_after_entrainment_forced=self.locals.cloud_total_water_after_entrainment_forced,
                    cloud_total_water_after_entrainment_downdraft_forced=self.locals.cloud_total_water_after_entrainment_downdraft_forced,
                    downdraft_saturation_vapor_forced=self.locals.downdraft_saturation_vapor_forced,
                    vapor_wetbulb=self.locals.vapor_wetbulb,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    environment_moist_static_energy_forced=self.locals.environment_moist_static_energy_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    cloud_moist_static_energy_downdraft_forced=self.locals.cloud_moist_static_energy_downdraft_forced,
                    evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                    geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                    mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    gamma_cloud_levels_forced=self.locals.gamma_cloud_levels_forced,
                    total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                    total_normalized_integrated_evaporate_forced=self.locals.total_normalized_integrated_evaporate_forced,
                    buoyancy=self.locals.buoyancy,
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
                    geopotential_height_cloud_levels=self.locals.geopotential_height_cloud_levels,
                    geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_updraft=self.locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    d_buoyancy=self.locals.d_buoyancy,
                    d_buoyancy_forced=self.locals.d_buoyancy_forced,
                    gamma_cloud_levels=self.locals.gamma_cloud_levels,
                    gamma_cloud_levels_forced=self.locals.gamma_cloud_levels_forced,
                    t_cloud_levels=self.locals.t_cloud_levels,
                    t_cloud_levels_forced=self.locals.t_cloud_levels_forced,
                    cloud_workfunction_0=self.locals.cloud_workfunction_0,
                    cloud_workfunction_1=self.locals.cloud_workfunction_1,
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
                    geopotential_height_cloud_levels=self.locals.geopotential_height_cloud_levels,
                    geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_updraft=self.locals.normalized_massflux_updraft,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    d_buoyancy=self.locals.d_buoyancy,
                    d_buoyancy_forced=self.locals.d_buoyancy_forced,
                    gamma_cloud_levels=self.locals.gamma_cloud_levels,
                    gamma_cloud_levels_forced=self.locals.gamma_cloud_levels_forced,
                    t_cloud_levels=self.locals.t_cloud_levels,
                    t_cloud_levels_forced=self.locals.t_cloud_levels_forced,
                    cin_0=self.locals.cin_0,
                    cin_1=self.locals.cin_1,
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
                    cin_0=self.locals.cin_0,
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
                    geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                    t_old=state.input_output.t_old,
                    t_new=self.locals.t_new,
                    t_cloud_levels_forced=self.locals.t_cloud_levels_forced,
                    vapor_old=state.input_output.vapor_old,
                    vapor_forced=self.locals.vapor_forced,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    vertical_velocity_2d=self.locals.vertical_velocity_2d,
                    cape_removal_time_scale=self.locals.cape_removal_time_scale,
                    cape_removal_time_scale_from_state=state.output.cape_removal_time_scale,
                    pbl_time_scale=self.locals.pbl_time_scale,
                    pbl_time_scale_from_state=state.output.pbl_time_scale,
                    cloud_work_function_1_pbl=self.locals.cloud_workfunction_1_pbl,
                    cloud_work_function_1_fa=self.locals.cloud_workfunction_1_fa,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # Trigger function based on Xie et al 2019
                # NOTE not implemented, does not run with test config
                self._Xie_trigger_function(plume_dependent_constants=self.plume_dependent_constants)

                # determine downdraft strength in terms of windshear
                # NOTE test GF2020_CumulusParameterization_DowndraftWindShear_{plume}:
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
                    psum=self.locals.psum,
                    psumh=self.locals.psumh,
                    total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                    total_normalized_integrated_evaporate_forced=self.locals.total_normalized_integrated_evaporate_forced,
                    scale_dependence_factor_downdraft=self.locals.scale_dependence_factor_downdraft,
                    epsilon=self.locals.epsilon,
                    epsilon_min=self.locals.epsilon_min,
                    epsilon_max=self.locals.epsilon_max,
                    epsilon_computed=self.locals.epsilon_computed,
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
                    environment_massflux=self.locals.environment_massflux,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # check mass conservation
                # NOTE This code runs in the Fortran and only has one output: totmas (total mass).
                # totmas has only one use: a conditional log write if total mass is above a 1e-6
                # with a disabled (commented) fatal error call.
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
                    geopotential_height_cloud_levels_forced=self.locals.geopotential_height_cloud_levels_forced,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    environment_massflux=self.locals.environment_massflux,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    c1d=self.locals.c1d,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    u_cloud_levels=self.locals.u_cloud_levels,
                    v_cloud_levels=self.locals.v_cloud_levels,
                    u_c=self.locals.u_c,
                    v_c=self.locals.v_c,
                    u_c_downdraft=self.locals.u_c_downdraft,
                    v_c_downdraft=self.locals.v_c_downdraft,
                    cloud_moist_static_energy_forced=self.locals.cloud_moist_static_energy_forced,
                    cloud_moist_static_energy_downdraft_forced=self.locals.cloud_moist_static_energy_downdraft_forced,
                    environment_moist_static_energy_cloud_levels_forced=self.locals.environment_moist_static_energy_cloud_levels_forced,
                    vapor_cloud_levels_forced=self.locals.vapor_cloud_levels_forced,
                    cloud_total_water_after_entrainment_forced=self.locals.cloud_total_water_after_entrainment_forced,
                    cloud_total_water_after_entrainment_downdraft_forced=self.locals.cloud_total_water_after_entrainment_downdraft_forced,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                    melting=self.locals.melting,
                    partition_liquid_ice=self.locals.partition_liquid_ice,
                    epsilon_forced=state.output.epsilon_forced,
                    d_buoyancy_downdraft_forced=self.locals.d_buoyancy_downdraft_forced,
                    del_u_cloud_ensemble=self.locals.del_u_cloud_ensemble,
                    del_v_cloud_ensemble=self.locals.del_v_cloud_ensemble,
                    del_moist_static_energy_cloud_ensemble=self.locals.del_moist_static_energy_cloud_ensemble,
                    del_t_cloud_ensemble=self.locals.del_t_cloud_ensemble,
                    del_vapor_cloud_ensemble=self.locals.del_vapor_cloud_ensemble,
                    del_cloud_liquid_cloud_ensemble=self.locals.del_cloud_liquid_cloud_ensemble,
                    del_buoyancy_cloud_ensemble=self.locals.del_buoyancy_cloud_ensemble,
                    t_tendency_from_environmental_subsidence=self.locals.t_tendency_from_environmental_subsidence,
                    moist_static_energy_tendency_from_environmental_subsidence=self.locals.moist_static_energy_tendency_from_environmental_subsidence,
                    vapor_tendency_from_environmental_subsidence=self.locals.vapor_tendency_from_environmental_subsidence,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # apply environmental subsidence on grid-scale ice and
                # liq water contents, and cloud fraction (Upwind scheme)
                # NOTE not implemented, does not run with test config
                self._environmental_subsidence()

                # make the smoothness procedure
                # NOTE test GF2020_CumulusParameterization_SmoothTendencies_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._smooth_tendencies(
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    del_moist_static_energy_cloud_ensemble=self.locals.del_moist_static_energy_cloud_ensemble,
                    del_vapor_cloud_ensemble=self.locals.del_vapor_cloud_ensemble,
                    del_cloud_liquid_cloud_ensemble=self.locals.del_cloud_liquid_cloud_ensemble,
                    del_u_cloud_ensemble=self.locals.del_u_cloud_ensemble,
                    del_v_cloud_ensemble=self.locals.del_v_cloud_ensemble,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # using smoothed tendencies, calculate changed environmental profiles
                # NOTE test GF2020_CumulusParameterization_ModifyEnvironmentProfiles_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._modify_environment_profiles(
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    updraft_origin_level=state.output.updraft_origin_level,
                    ocean_fraction=state.input.ocean_fraction,
                    p_forced=state.input_output.p_forced,
                    t_new=self.locals.t_new,
                    t_modified=self.locals.t_modified,
                    vapor_forced=self.locals.vapor_forced,
                    vapor_modified=self.locals.vapor_modified,
                    environment_moist_static_energy_forced=self.locals.environment_moist_static_energy_forced,
                    environment_moist_static_energy_modified=self.locals.environment_moist_static_energy_modified,
                    moist_static_energy_origin_level_forced=self.locals.moist_static_energy_origin_level_forced,
                    moist_static_energy_origin_level_modified=self.locals.moist_static_energy_origin_level_modified,
                    partition_liquid_ice=self.locals.partition_liquid_ice,
                    del_moist_static_energy_cloud_ensemble=self.locals.del_moist_static_energy_cloud_ensemble,
                    del_t_cloud_ensemble=self.locals.del_t_cloud_ensemble,
                    del_vapor_cloud_ensemble=self.locals.del_vapor_cloud_ensemble,
                    del_cloud_liquid_cloud_ensemble=self.locals.del_cloud_liquid_cloud_ensemble,
                    del_u_cloud_ensemble=self.locals.del_u_cloud_ensemble,
                    del_v_cloud_ensemble=self.locals.del_v_cloud_ensemble,
                    moist_static_energy_tendency_from_environmental_subsidence=self.locals.moist_static_energy_tendency_from_environmental_subsidence,
                    vapor_tendency_from_environmental_subsidence=self.locals.vapor_tendency_from_environmental_subsidence,
                    t_tendency_from_environmental_subsidence=self.locals.t_tendency_from_environmental_subsidence,
                    arbitrary_numerical_parameter=self.locals.arbitrary_numerical_parameter,
                    AVERAGE_LAYER_DEPTH=self.plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # calculate moist static energy, heights, environmental saturation mixing ratio
                # NOTE test GF2020_CumulusParameterization_EnvironmentConditions_3_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ❌ one field, one point (0.0%), 2 ULP
                # NOTE      shallow ✅
                self._environment_conditions(
                    p=state.input_output.p_forced,
                    p_surface=state.input_output.p_surface,
                    t=self.locals.t_modified,
                    vapor=self.locals.vapor_modified,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    moist_static_energy=self.locals.environment_moist_static_energy_modified,
                    saturation_moist_static_energy=self.locals.environment_saturation_moist_static_energy_modified,
                    saturation_mixing_ratio=self.locals.environment_saturation_mixing_ratio_modified,
                    geopotential_height=self.locals.geopotential_height_modified,
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # environmental values on cloud levels
                # NOTE test GF2020_CumulusParameterization_EnvironmentCloudLevels_3_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ❌ one field, one point (0.0%), 32 ULP
                # NOTE      shallow ❌ one field, one point (0.0%), 32 ULP
                p_3d = state.output.p_cloud_levels_forced.field[
                    :, :, :, self.plume_dependent_constants.PLUME_INDEX
                ]
                self._environment_cloud_levels(
                    p=state.input_output.p_forced,
                    p_surface=state.input_output.p_surface,
                    p_cloud_levels=p_3d,
                    topography_height_no_negative=state.input_output.topography_height_no_negative,
                    geopotential_height=self.locals.geopotential_height_modified,
                    geopotential_height_cloud_levels=self.locals.geopotential_height_cloud_levels_modified,
                    t=self.locals.t_modified,
                    t_surface=state.input_output.t_surface,
                    t_cloud_levels=self.locals.t_cloud_levels_modified,
                    vapor=self.locals.vapor_modified,
                    vapor_cloud_levels=self.locals.vapor_cloud_levels_modified,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    u_cloud_levels=self.locals.u_cloud_levels,
                    v_cloud_levels=self.locals.v_cloud_levels,
                    environment_saturation_mixing_ratio=self.locals.environment_saturation_mixing_ratio_modified,
                    environment_saturation_mixing_ratio_cloud_levels=self.locals.environment_saturation_mixing_ratio_cloud_levels_modified,
                    environment_moist_static_energy=self.locals.environment_moist_static_energy_modified,
                    environment_moist_static_energy_cloud_levels=self.locals.environment_moist_static_energy_cloud_levels_modified,
                    environment_saturation_moist_static_energy=self.locals.environment_saturation_moist_static_energy_modified,
                    environment_saturation_moist_static_energy_cloud_levels=self.locals.environment_saturation_moist_static_energy_cloud_levels_modified,
                    gamma_cloud_levels=self.locals.gamma_cloud_levels,
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )
                state.output.p_cloud_levels_forced.field[
                    :, :, :, self.plume_dependent_constants.PLUME_INDEX
                ] = p_3d

                # static control
                # NOTE test GF2020_CumulusParameterization_StaticControl_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._static_control(
                    error_code=state.output.error_code,
                    start_level=self.locals.start_level,
                    lcl_level=state.output.lcl_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    cloud_moist_static_energy_modified=self.locals.cloud_moist_static_energy_modified,
                    moist_static_energy_origin_level_modified=self.locals.moist_static_energy_origin_level_modified,
                    environment_saturation_moist_static_energy_modified=self.locals.environment_saturation_moist_static_energy_modified,
                    environment_moist_static_energy_cloud_levels_modified=self.locals.environment_moist_static_energy_cloud_levels_modified,
                    environment_saturation_moist_static_energy_cloud_levels_modified=self.locals.environment_saturation_moist_static_energy_cloud_levels_modified,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                    normalized_massflux_updraft_modified=self.locals.normalized_massflux_updraft_modified,
                    partition_liquid_ice=self.locals.partition_liquid_ice,
                    vapor_excess=self.locals.vapor_excess,
                    t_excess=self.locals.t_excess,
                    add_buoyancy=self.locals.add_buoyancy,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    d_buoyancy_modified=self.locals.d_buoyancy_modified,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # workfunctions for updraft
                # NOTE test GF2020_CumulusParameterization_UpdraftWorkfunctions_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._updraft_workfunctions(
                    error_code=state.output.error_code,
                    updraft_origin_level=state.output.updraft_origin_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    geopotential_height_cloud_levels_modified=self.locals.geopotential_height_cloud_levels_modified,
                    normalized_massflux_updraft_modified=self.locals.normalized_massflux_updraft_modified,
                    d_buoyancy_modified=self.locals.d_buoyancy_modified,
                    gamma_cloud_levels=self.locals.gamma_cloud_levels,
                    t_cloud_levels_modified=self.locals.t_cloud_levels_modified,
                    cloud_workfunction_0_modified=self.locals.cloud_workfunction_0_modified,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                    epsilon_forced=state.output.epsilon_forced,
                    precipitation_ensemble=self.locals.precipitation_ensemble,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # large scale forcing
                # calculate cloud base mass flux
                # NOTE test GF2020_CumulusParameterization_CloudBaseMassFlux_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._large_scale_forcing(
                    error_code=state.output.error_code,
                    error_code_2=self.locals.error_code_2,
                    error_code_3=self.locals.error_code_3,
                    updraft_origin_level=state.output.updraft_origin_level,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    pbl_level=state.input_output.pbl_level,
                    ocean_fraction=self.locals.ocean_fraction,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    vapor_forced=self.locals.vapor_forced,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    effective_condensate_to_fall_forced=self.locals.effective_condensate_to_fall_forced,
                    evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                    omega=state.input_output.omega,
                    convective_scale_velocity=state.input_output.convective_scale_velocity,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    cloud_moist_static_energy=self.locals.cloud_moist_static_energy,
                    cloud_moist_static_energy_forced=self.locals.cloud_moist_static_energy_forced,
                    environment_moist_static_energy_cloud_levels=self.locals.environment_moist_static_energy_cloud_levels,
                    environment_moist_static_energy_cloud_levels_forced=self.locals.environment_moist_static_energy_cloud_levels_forced,
                    dmoist_static_energydt=self.locals.dmoist_static_energydt,
                    cloud_workfunction_0=self.locals.cloud_workfunction_0,
                    cloud_workfunction_0_modified=self.locals.cloud_workfunction_0_modified,
                    cloud_workfunction_1=self.locals.cloud_workfunction_1,
                    cloud_workfunction_1_pbl=self.locals.cloud_workfunction_1_pbl,
                    arbitrary_numerical_parameter=self.locals.arbitrary_numerical_parameter,
                    f_dicycle_modified=self.locals.f_dicycle_modified,
                    cape_removal_time_scale=self.locals.cape_removal_time_scale,
                    epsilon_forced=state.output.epsilon_forced,
                    k_x_modified=self.locals.k_x_modified,
                    mass_flux_ensemble=self.locals.mass_flux_ensemble,
                    precipitation_ensemble=self.locals.precipitation_ensemble,
                    xff_mid=self.locals.xff_mid,
                    plume_dependent_constants=self.plume_dependent_constants,
                )

                # Include kinetic energy dissipation converted to heating
                # NOTE test GF2020_CumulusParameterization_KeToHeating_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._kinetic_energy_to_heating(
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    u=state.input_output.u,
                    v=state.input_output.v,
                    del_u_cloud_ensemble=self.locals.del_u_cloud_ensemble,
                    del_v_cloud_ensemble=self.locals.del_v_cloud_ensemble,
                    del_t_cloud_ensemble=self.locals.del_t_cloud_ensemble,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # feedback
                self._feedback()

                # net precipitation flux (after downdraft evaporation)
                # NOTE test GF2020_CumulusParameterization_PrecipitationFlux_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._precipitation_flux(
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    cloud_base_mass_flux_modified=state.output.cloud_base_mass_flux_modified,
                    epsilon_forced=state.output.epsilon_forced,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                    precipitation_flux=self.locals.precipitation_flux,
                    evaporation_flux=self.locals.evaporation_flux,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # rainfall evap below cloud base
                # NOTE test GF2020_CumulusParameterization_RainEvaporationBelowCloudBase_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._rain_evaporation_below_cloud_base(
                    error_code=state.output.error_code,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    ocean_fraction=state.input.ocean_fraction,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    p_surface=state.input_output.p_surface,
                    t_cloud_levels=self.locals.t_cloud_levels,
                    vapor_cloud_levels_forced=self.locals.vapor_cloud_levels_forced,
                    environment_saturation_mixing_ratio_cloud_levels=self.locals.environment_saturation_mixing_ratio_cloud_levels,
                    epsilon_forced=state.output.epsilon_forced,
                    cloud_base_mass_flux_modified=state.output.cloud_base_mass_flux_modified,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                    precip=state.output.precip,
                    precipitation_flux=self.locals.precipitation_flux,
                    evaporation_flux=self.locals.evaporation_flux,
                    evaporation_below_cloud_base=self.locals.evaporation_below_cloud_base,
                    dtdt=state.output.dtdt,
                    dvapordt=state.output.dvapordt,
                    dbuoyancydt=state.output.dbuoyancydt,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # includes effects of the remained cloud dissipation into the enviroment
                # NOTE test GF2020_CumulusParameterization_CloudDissipation_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._cloud_dissapation(
                    error_code=state.output.error_code,
                    updraft_lfc_level=state.output.updraft_lfc_level,
                    cloud_top_level=state.output.cloud_top_level,
                    hydrostatic_air_density=self.locals.hydrostatic_air_density,
                    geopotential_height_forced=state.input_output.geopotential_height_forced,
                    t_cloud_levels_forced=self.locals.t_cloud_levels_forced,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    cloud_base_mass_flux_modified=state.output.cloud_base_mass_flux_modified,
                    vapor_cloud_levels_forced=self.locals.vapor_cloud_levels_forced,
                    cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                    environment_saturation_mixing_ratio_cloud_levels_forced=self.locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                    environment_saturation_moist_static_energy_cloud_levels_forced=self.locals.environment_saturation_moist_static_energy_cloud_levels_forced,
                    vertical_velocity_3d=self.locals.vertical_velocity_3d,
                    scale_dependence_factor=state.output.scale_dependence_factor,
                    dtdt=state.output.dtdt,
                    dvapordt=state.output.dvapordt,
                    dcloudicedt=state.output.dcloudicedt,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # total (deep+mid) evaporation flux for output (units kg/kg/s)
                # NOTE test GF2020_CumulusParameterization_TotalEvaporationFlux_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._total_evaporation_flux(
                    error_code=state.output.error_code,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                    cloud_top_level=state.output.cloud_top_level,
                    p_cloud_levels_forced=state.output.p_cloud_levels_forced,
                    evaporation_flux=self.locals.evaporation_flux,
                    evaporation_sublimation_tendency=state.output.evaporation_sublimation_tendency,
                )

                # lightning flashes density (parameterization from Lopez 2016, MWR)
                # NOTE this section does not run in the test case, and has not been implemented.
                self._lightning_flash_density()

                # output precipitation (only deep plume)
                # NOTE test GF2020_CumulusParameterization_DeepPrecipitationOutput_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._deep_precipitation_output(
                    error_code=state.output.error_code,
                    cloud_top_level=state.output.cloud_top_level,
                    precipitation_flux=self.locals.precipitation_flux,
                    convective_precip_flux=state.output.convective_precip_flux,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # for tracer convective transport / outputs
                # NOTE test GF2020_CumulusParameterization_TracerOutput_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._tracer_output(
                    error_code=state.output.error_code,
                    updraft_column_temperature_forced=self.locals.updraft_column_temperature_forced,
                    t_cloud_levels=self.locals.t_cloud_levels,
                    t_updraft=state.output.t_updraft,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # convert mass fluxes, etc...
                # NOTE test GF2020_CumulusParameterization_PrepareOutput_{plume}:
                # NOTE      deep ✅
                # NOTE      mid ✅
                # NOTE      shallow ✅
                self._prepare_output(
                    error_code=state.output.error_code,
                    cloud_base_mass_flux_modified=state.output.cloud_base_mass_flux_modified,
                    total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                    total_normalized_integrated_evaporate_forced=self.locals.total_normalized_integrated_evaporate_forced,
                    normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                    normalized_massflux_downdraft_forced=state.output.normalized_massflux_downdraft_forced,
                    condensate_to_fall_forced=state.output.condensate_to_fall_forced,
                    evaporate_in_downdraft_forced=state.output.evaporate_in_downdraft_forced,
                    mass_entrainment_updraft_forced=state.output.mass_entrainment_updraft_forced,
                    mass_detrainment_updraft_forced=state.output.mass_detrainment_updraft_forced,
                    mass_entrainment_downdraft_forced=state.output.mass_entrainment_downdraft_forced,
                    mass_detrainment_downdraft_forced=state.output.mass_detrainment_downdraft_forced,
                    environment_massflux=self.locals.environment_massflux,
                    vapor_tendency_from_environmental_subsidence=self.locals.vapor_tendency_from_environmental_subsidence,
                    moist_static_energy_tendency_from_environmental_subsidence=self.locals.moist_static_energy_tendency_from_environmental_subsidence,
                    t_tendency_from_environmental_subsidence=self.locals.t_tendency_from_environmental_subsidence,
                    plume=self.plume_dependent_constants.PLUME_INDEX,
                )

                # outputs a model sounding for the stand-alone code (part 2)
                # NOTE this section does not run in the test case, and has not been implemented.
                self._sounding()

                # section for atmospheric composition
                # NOTE this section does not run in the test case, and has not been implemented.
                self._atmospheric_composition()

                # begin: for GATE soundings
                # NOTE this section does not run in the test case, and has not been implemented.
                self._gate_sounding()
