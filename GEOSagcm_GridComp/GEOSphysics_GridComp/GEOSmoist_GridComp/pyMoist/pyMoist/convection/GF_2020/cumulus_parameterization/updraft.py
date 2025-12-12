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
from ndsl.dsl.typing import (
    FloatField,
    FloatFieldIJ,
    IntFieldIJ,
    Int,
)
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import IntFieldIJ_Plume


def updraft_air_temperature(
    error_code: IntFieldIJ_Plume,
    updraft_t: FloatField,
    cloud_moist_static_energy_forced: FloatField,
    geopotential_height_cloud_levels_forced: FloatField,
    cloud_vapor_mixing_ratio_forced: FloatField,
    t_cloud_levels_forced: FloatField,
    plume: Int,
):
    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            updraft_t = (1.0 / cumulus_parameterization_constants.CP) * (
                cloud_moist_static_energy_forced
                - constants.MAPL_GRAV * geopotential_height_cloud_levels_forced
                - cumulus_parameterization_constants.XLV * cloud_vapor_mixing_ratio_forced
            )

    with computation(PARALLEL), interval(-1, None):
        if error_code[0, 0][plume] == 0:
            updraft_t = t_cloud_levels_forced

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] != 0:
            updraft_t = t_cloud_levels_forced


def cup_up_aa0(
    local_buoyancy: FloatField,
    local_gamma_cloud_levels: FloatField,
    cloud_top_level: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    updraft_origin_level: IntFieldIJ_Plume,
    local_geopotential_height_cloud_levels: FloatField,
    local_t_cloud_levels: FloatField,
    local_normalized_massflux_updraft: FloatField,
    local_integ: IntFieldIJ,
    local_integ_interval: IntFieldIJ,
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_cloud_work_function: FloatFieldIJ,
):
    from __externals__ import k_start

    with computation(FORWARD), interval(...):
        local_cloud_work_function = 0.0

        if local_integ == 1:
            if local_integ_interval == cumulus_parameterization_constants.BL:
                kbeg = k_start
                kend = updraft_lfc_level[0, 0][plume] - 2
            elif local_integ_interval == cumulus_parameterization_constants.CIN:
                kbeg = updraft_origin_level[0, 0][plume] - 1
                kend = updraft_lfc_level[0, 0][plume] - 2

        else:
            kbeg = updraft_lfc_level[0, 0][plume] - 1
            kend = cloud_top_level[0, 0][plume] - 1

    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if K >= kbeg and K <= kend:
                dz = local_geopotential_height_cloud_levels[0, 0, 1] - local_geopotential_height_cloud_levels
                aa_1 = (
                    local_normalized_massflux_updraft
                    * (constants.MAPL_GRAV / (cumulus_parameterization_constants.CP * local_t_cloud_levels))
                    * local_buoyancy
                    / (1.0 + local_gamma_cloud_levels)
                )
                aa_2 = (
                    local_normalized_massflux_updraft[0, 0, 1]
                    * (
                        constants.MAPL_GRAV
                        / (cumulus_parameterization_constants.CP * local_t_cloud_levels[0, 0, 1])
                    )
                    * local_buoyancy[0, 0, 1]
                    / (1.0 + local_gamma_cloud_levels[0, 0, 1])
                )
                da = 0.5 * (aa_1 + aa_2) * dz

                local_cloud_work_function = local_cloud_work_function + da


def cloud_work_function_zero(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_cloud_work_function: FloatFieldIJ,
):
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if local_cloud_work_function == 0.0:
                error_code[0, 0][plume] = 17
                # ierrc[0,0][plume]="cloud work function zero"


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
        #     cloud_top_level=state.output.cloud_top_level,
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
            cloud_top_level=state.output.cloud_top_level,
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
            cloud_top_level=state.output.cloud_top_level,
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
        self._cup_up_aa0(
            local_buoyancy=locals.buoyancy,
            local_gamma_cloud_levels=locals.gamma_cloud_levels,
            cloud_top_level=state.output.cloud_top_level,
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
            cloud_top_level=state.output.cloud_top_level,
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


class UpdraftUpdateWorkfunctions:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
