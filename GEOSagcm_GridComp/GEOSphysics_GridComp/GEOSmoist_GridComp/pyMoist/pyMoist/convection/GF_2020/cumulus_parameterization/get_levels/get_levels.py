from ndsl import StencilFactory, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.get_levels.stencils import (
    find_maximum_updraft_origin_level,
    find_detrainmet_start_level,
    find_highest_moist_static_energy_level,
    find_lcl,
    set_start_level,
    convective_cloud_base_level,
    updraft_rates_pdf,
    cloud_top_checks,
)
from ndsl.logging import ndsl_log


class MaximumUpdraftOriginLevel:
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
        self._find_maximum_updraft_origin_level = stencil_factory.from_dims_halo(
            func=find_maximum_updraft_origin_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._find_maximum_updraft_origin_level(
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            topography_height_no_negative=state.input_output.topography_height_no_negative,
            error_code=state.output.error_code,
            maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
            MAX_UPDRAFT_ORIGIN_HEIGHT=plume_dependent_constants.MAX_UPDRAFT_ORIGIN_HEIGHT,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class DowndraftDetrainmentLevel:
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
        self._find_detrainmet_start_level = stencil_factory.from_dims_halo(
            func=find_detrainmet_start_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._find_detrainmet_start_level(
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            topography_height_no_negative=state.input_output.topography_height_no_negative,
            error_code=state.output.error_code,
            detrainment_start_level=locals.detrainment_start_level,
            DETRAINMENT_CRITICAL_DEPTH=plume_dependent_constants.DETRAINMENT_CRITICAL_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class HighestMoistStaticEnergyLevel:
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
        self._find_highest_moist_static_energy_level = stencil_factory.from_dims_halo(
            func=find_highest_moist_static_energy_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._find_highest_moist_static_energy_level(
            moist_static_energy=locals.environment_moist_static_energy_cloud_levels_forced,
            error_code=state.output.error_code,
            maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
            updraft_origin_level=state.output.updraft_origin_level,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class GetLCL:
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
        self._find_lcl = stencil_factory.from_dims_halo(
            func=find_lcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
                "ADV_TRIGGER": config.ADV_TRIGGER,
            },
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        if self.cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called GetLCL with "
                "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
            )

        if self.config.ADV_TRIGGER != 1:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called GetLCL with "
                "untested ADV_TRIGGER option. Running untested code... proceed with caution"
            )

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
            AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class ConvectiveCloudBaseLevel:
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
        self._set_start_level = stencil_factory.from_dims_halo(
            func=set_start_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        self._convective_cloud_base_level = stencil_factory.from_dims_halo(
            func=convective_cloud_base_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "OVERSHOOT": cumulus_parameterization_config.OVERSHOOT,
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                "MOIST_TRIGGER": cumulus_parameterization_config.MOIST_TRIGGER,
                "USE_MEMORY": cumulus_parameterization_config.USE_MEMORY,
                "BOUNDARY_CONDITION_METHOD": cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD,
            },
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        if self.cumulus_parameterization_config.OVERSHOOT != 0:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called ConvectiveCloudBaseLevel with "
                "untested ZERO_DIFF option. Running untested code... proceed with caution"
            )

        if self.cumulus_parameterization_config.ZERO_DIFF != 0:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called ConvectiveCloudBaseLevel with "
                "untested ZERO_DIFF option. Running untested code... proceed with caution"
            )

        if self.cumulus_parameterization_config.MOIST_TRIGGER != 0:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called ConvectiveCloudBaseLevel with "
                "untested MOIST_TRIGGER option. Running untested code... proceed with caution"
            )

        if self.cumulus_parameterization_config.USE_MEMORY != -1:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called ConvectiveCloudBaseLevel with "
                "untested USE_MEMORY option. Running untested code... proceed with caution"
            )

        if self.cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called ConvectiveCloudBaseLevel with "
                "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
            )

        self._set_start_level(
            updraft_origin_level=state.output.updraft_origin_level,
            start_level=locals.start_level,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._convective_cloud_base_level(
            error_code=state.output.error_code,
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
            cloud_top=state.output.cloud_top,
            AVERAGE_LAYER_DEPTH=plume_dependent_constants.AVERAGE_LAYER_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class CloudTop:
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
        self._updraft_rates_pdf = stencil_factory.from_dims_halo(
            func=updraft_rates_pdf,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"OVERSHOOT": cumulus_parameterization_config.OVERSHOOT},
        )

        self._cloud_top_checks = stencil_factory.from_dims_halo(
            func=cloud_top_checks,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        if self.cumulus_parameterization_config.OVERSHOOT != 0:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called CloudTop with "
                "untested OVERSHOOT option. Running untested code... proceed with caution"
            )

        self._updraft_rates_pdf(
            entrainment_rate=state.output.entrainment_rate,
            moist_static_energy=locals.environment_moist_static_energy_forced,
            saturation_moist_static_energy=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
            moist_static_energy_origin_level=locals.moist_static_energy_origin_level_forced,
            updraft_lfc_level=state.output.updraft_lfc_level,
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            cloud_moist_static_energy=locals.cloud_moist_static_energy_forced_transported,
            error_code=state.output.error_code,
            cloud_top=state.output.cloud_top,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._cloud_top_checks(
            cloud_top=state.output.cloud_top,
            p=state.output.p_cloud_levels_forced,
            geopotential_height=locals.geopotential_height_cloud_levels,
            error_code=state.output.error_code,
            last_error_code=state.input.last_error_code,
            updraft_lfc_level=state.output.updraft_lfc_level,
            MINIMUM_DEPTH=plume_dependent_constants.MINIMUM_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class DowndraftOriginLevel:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
