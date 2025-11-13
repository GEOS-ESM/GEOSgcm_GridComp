from ndsl import StencilFactory, QuantityFactory, ndsl_log
from pyMoist.convection.GF_2020.cumulus_parameterization.environment.stencils import (
    environment_conditions,
    environment_cloud_levels,
    environment_mass_flux,
)
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import (
    GF2020CumulusParameterizationConfig,
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.cumulus_parameterization.state import (
    GF2020CumulusParameterizationState,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import (
    GF2020CumulusParameterizationLocals,
)
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)


class EnvironmentConditions:
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
        self._environment_conditions = stencil_factory.from_dims_halo(
            func=environment_conditions,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "SATURATION_CALCULATION_CHOICE": cumulus_parameterization_config.SATURATION_CALCULATION_CHOICE
            },
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
        data_type: int,
    ):
        if self.cumulus_parameterization_config.SATURATION_CALCULATION_CHOICE != 1:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called EnvironmentalConditions with "
                "untested SATURATION_CALCULATION_CHOICE option. Running untested code... proceed with caution"
            )

        if data_type == 0:
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
                plume=plume_dependent_constants.PLUME_INDEX,
            )
        elif data_type == 1:
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
                plume=plume_dependent_constants.PLUME_INDEX,
            )
        else:
            raise NotImplementedError("EnvironmentCloudLevels call type not supported.")


class EnvironmentCloudLevels:
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
        self._environment_cloud_levels = stencil_factory.from_dims_halo(
            func=environment_cloud_levels,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "CLOUD_LEVEL_GRID": cumulus_parameterization_config.CLOUD_LEVEL_GRID
            },
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
        data_type: int,
    ):

        if self.cumulus_parameterization_config.CLOUD_LEVEL_GRID != 1:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called EnvironmentalCloudLevels with "
                "untested CLOUD_LEVEL_GRID option. Running untested code... proceed with caution"
            )

        if data_type == 0:
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
                plume=plume_dependent_constants.PLUME_INDEX,
            )
        elif data_type == 1:
            p_3d = state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ]
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
                vapor_cloud_levels=locals.vapor_cloud_levels,
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
                plume=plume_dependent_constants.PLUME_INDEX,
            )
            state.output.p_cloud_levels_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ] = p_3d
        else:
            raise NotImplementedError("EnvironmentCloudLevels call type not supported.")


class EnvironmentMassFlux:
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
        self._environment_mass_flux = stencil_factory.from_dims_halo(
            func=environment_mass_flux,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._environment_mass_flux(
            # zenv=,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            # zuo=,
            # edto=,
            # zdo=,
        )


class EnvironmentalSubsidence:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
