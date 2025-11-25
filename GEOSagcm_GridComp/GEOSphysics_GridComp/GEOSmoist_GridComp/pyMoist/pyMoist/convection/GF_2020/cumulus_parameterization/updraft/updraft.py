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
from pyMoist.convection.GF_2020.cumulus_parameterization.updraft.stencils import (
    cup_up_aa0,
    cloud_work_function_zero,
    in_cloud_updraft_air_temperature,
)


class UpdraftMassFluxProfile:
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
        # self._updraft_rates_pdf = stencil_factory.from_dims_halo(
        #     func=updraft_rates_pdf,
        #     compute_dims=[X_DIM, Y_DIM, Z_DIM],
        #     externals={"OVERSHOOT": cumulus_parameterization_config.OVERSHOOT},
        # )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        pass
        # self._updraft_rates_pdf(
        #     entrainment_rate=state.output.entrainment_rate,
        #     moist_static_energy=locals.environment_moist_static_energy_forced,
        #     saturation_moist_static_energy=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
        #     moist_static_energy_origin_level=locals.moist_static_energy_origin_level_forced,
        #     updraft_lfc_level=state.output.updraft_lfc_level,
        #     geopotential_height=locals.geopotential_height_cloud_levels_forced,
        #     cloud_moist_static_energy=locals.cloud_moist_static_energy_forced_transported,
        #     error_code=state.output.error_code,
        #     cloud_top=state.output.cloud_top,
        #     plume=plume_dependent_constants.PLUME_INDEX,
        # )


class UpdraftMoistureProfile:
    def __init__():
        pass

    def __call__(self, *args, **kwds):
        pass


class UpdraftMoistStaticEnergyAndMomentumBudget:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class UpdraftInCloudUpdraftAirTemperature:
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
            externals={"FIRST_GUESS_W": cumulus_parameterization_config.FIRST_GUESS_W},
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._in_cloud_updraft_air_temperature(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_incloud_air_temp=locals.incloud_air_temp,
            local_cloud_moist_static_energy_forced=locals.cloud_moist_static_energy_forced,
            local_geopotential_height_cloud_levels_forced=locals.geopotential_height_cloud_levels_forced,
            local_incloud_water_vapor_mixing_ratio=locals.incloud_water_vapor_mixing_ratio,
            local_t_cloud_levels_forced=locals.t_cloud_levels_forced,
        )


class UpdraftInitialWorkfunctions:
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
        self._cup_up_aa0 = stencil_factory.from_dims_halo(
            func=cup_up_aa0,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._cloud_work_function_zero = stencil_factory.from_dims_halo(
            func=cloud_work_function_zero,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._cup_up_aa0(
            local_buoyancy=locals.buoyancy,
            local_gamma_cloud_levels=locals.gamma_cloud_levels,
            cloud_top=state.output.cloud_top,
            updraft_lfc_level=state.output.updraft_lfc_level,
            updraft_origin_level=state.output.updraft_origin_level,
            local_geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
            local_t_cloud_levels=locals.t_cloud_levels,
            local_normalized_massflux_updraft=locals.normalized_massflux_updraft,
            local_integ=locals.integ,
            local_integ_interval=locals.integ_interval,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_cloud_work_function=locals.cloud_work_function_0,
        )

        self._cup_up_aa0(
            local_buoyancy=locals.buoyancy_forced,
            local_gamma_cloud_levels=locals.gamma_cloud_levels_forced,
            cloud_top=state.output.cloud_top,
            updraft_lfc_level=state.output.updraft_lfc_level,
            updraft_origin_level=state.output.updraft_origin_level,
            local_geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels_forced,
            local_t_cloud_levels=locals.t_cloud_levels_forced,
            local_normalized_massflux_updraft=locals.normalized_massflux_updraft_forced,
            local_integ=locals.integ,
            local_integ_interval=locals.integ_interval,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_cloud_work_function=locals.cloud_work_function_1,
        )

        self._cloud_work_function_zero(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_cloud_work_function=locals.cloud_work_function_1,
        )


class UpdraftCIN:
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
        self._cup_up_aa0 = stencil_factory.from_dims_halo(
            func=cup_up_aa0,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        pass
        # self._cup_up_aa0(
        #     dby=dby,
        #     gamma_cup=gamma_cup,
        #     ktop=ktop,
        #     kbcon=kbcon,
        #     k22=k22,
        #     z_cup=z_cup,
        #     t_cup=t_cup,
        #     zu=zu,
        #     integ=integ,
        #     integ_interval=integ_interval,
        #     error_codd=error_code,
        #     plume=plume,
        #     aa0=cin0,
        # )

        # self._cup_up_aa0(
        #     dby=dbyo,
        #     gamma_cup=gammao_cup,
        #     ktop=ktop,
        #     kbcon=kbcon,
        #     k22=k22,
        #     z_cup=zo_cup,
        #     t_cup=tn_cup,
        #     zu=zuo,
        #     integ=integ,
        #     integ_interval=integ_interval,
        #     error_codd=error_code,
        #     plume=plume,
        #     aa0=cin1,
        # )


class UpdraftUpdateWorkfunctions:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
