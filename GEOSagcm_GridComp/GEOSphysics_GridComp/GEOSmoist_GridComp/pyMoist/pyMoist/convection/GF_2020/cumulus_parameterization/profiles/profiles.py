from ndsl import StencilFactory, QuantityFactory
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
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.profiles.stencils import (
    in_cloud_updraft_air_temperature,
    get_melting_profile,
    updraft_moisture_light,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import updraft_vertical_velocity


class C1DProfile:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._updraft_moisture_light = stencil_factory.from_dims_halo(
            func=updraft_moisture_light,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
                "MELT_GLAC": cumulus_parameterization_config.MELT_GLAC,
                "FRAC_MODIS": cumulus_parameterization_config.FRAC_MODIS,
                "QRC_CRIT_OCN": cumulus_parameterization_config.QRC_CRIT_OCN,
                "QRC_CRIT_LND": cumulus_parameterization_config.QRC_CRIT_LND,
            },
        )

        self._updraft_vertical_velocity = stencil_factory.from_dims_halo(
            func=updraft_vertical_velocity,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF},
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        if plume_dependent_constants.PLUME_INDEX == 2 and self.cumulus_parameterization_config.USE_C1D:
            # NOTE this needs to go in a stencil once NDSL base functions are merged
            locals.c1d.field[:] = self.config.C1

        if self.cumulus_parameterization_config.FIRST_GUESS_W or self.config.AUTOCONV == 4:
            self._updraft_moisture_light(
                start_level=locals.start_level,
                error_code=state.output.error_code,
                geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                cloud_vapor_mixing_ratio_forced=locals.cloud_vapor_mixing_ratio_forced,
                cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                precipitable_water_updraft_forced=state.output.precipitable_water_updraft_forced,
                total_normalized_integrated_condensate_forced=state.output.total_normalized_integrated_condensate_forced,
                cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
                unspecifid_temperature=locals.unspecifid_temperature,
                ocean_fraction=locals.ocean_fraction,
                convection_fraction=state.input.convection_fraction,
                surface_type=state.input.surface_type,
                p_forced=state.input_output.p_forced,
                cloud_top_level=state.output.cloud_top_level,
                buoyancy_forced=locals.buoyancy_forced,
                cloud_liquid_before_rain_forced=locals.cloud_liquid_before_rain_forced,
                t_cloud_levels=locals.t_cloud_levels,
                vapor_forced=locals.vapor_forced,
                gamma_cloud_levels_forced=locals.gamma_cloud_levels_forced,
                normalized_massflux_updraft_forced=state.output.normalized_massflux_updraft_forced,
                environment_saturation_mixing_ratio_cloud_levels_forced=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
                updraft_origin_level=state.output.updraft_origin_level,
                vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
                vapor_excess=locals.vapor_excess,
                mass_entrainment_updraft=locals.mass_entrainment_updraft,
                mass_detrainment_updraft=locals.mass_detrainment_updraft,
                psum=locals.psum,
                psumh=locals.psumh,
                c1d=locals.c1d,
                add_buoyancy=locals.add_buoyancy,
                AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
                C0=plume_dependent_constants.C0,
                plume=plume_dependent_constants.PLUME_INDEX,
            )

            self._updraft_vertical_velocity(
                vertical_velocity_3d=locals.vertical_velocity_3d,
                vertical_velocity_2d=locals.vertical_velocity_2d,
                convective_scale_velocity=state.input_output.convective_scale_velocity,
                entrainment_rate=state.output.entrainment_rate,
                detrainment_function_updraft=locals.detrainment_function_updraft,
                geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
                t_cloud_levels_forced=locals.t_cloud_levels_forced,
                unspecifid_temperature=locals.unspecifid_temperature,
                cloud_vapor_mixing_ratio_forced=locals.cloud_vapor_mixing_ratio_forced,
                cloud_liquid_after_rain_forced=state.output.cloud_liquid_after_rain_forced,
                vapor_forced=locals.vapor_forced,
                updraft_lfc_level=state.output.updraft_lfc_level,
                cloud_top_level=state.output.cloud_top_level,
                error_code=state.output.error_code,
                plume=plume_dependent_constants.PLUME_INDEX,
            )


class MeltingProfile:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._get_melting_profile = stencil_factory.from_dims_halo(
            func=get_melting_profile,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "MELT_GLAC": cumulus_parameterization_config.MELT_GLAC,
            },
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        self._get_melting_profile(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_melting_layer=locals.melting_layer,
            local_partition_liquid_ice=locals.partition_liquid_ice,
            p_cloud_levels_forced=state.output.p_cloud_levels_forced,
            precipitable_water_updraft_forced=state.output.precipitable_water_updraft_forced,
            local_melting=locals.melting,
        )


class InCloudTemperature:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._in_cloud_updraft_air_temperature = stencil_factory.from_dims_halo(
            func=in_cloud_updraft_air_temperature,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        pass
        self._in_cloud_updraft_air_temperature(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_incloud_air_temp_forced=locals.incloud_air_temp_forced,
            local_hcdo=locals.hcdo,
            local_geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
            local_incloud_water_vapor_mixing_ratio_forced=locals.cloud_vapor_mixing_ratio_forced,
            local_t_cloud_levels_forced=locals.t_cloud_levels_forced,
        )
